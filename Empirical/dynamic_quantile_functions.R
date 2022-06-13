Plotnetivrq <- function (x, parm = NULL, level = 0.95, mfrow = NULL, 
          mar = NULL, ylim = NULL, main = NULL, col = gray(c(0, 0.75)), 
          border = NULL, lcol = 2, lty = 1:2, cex = 0.5, pch = 20, 
          type = "b", xlab = "", ylab = "") 
{
  zalpha <- qnorm(1 - (1 - level)/2)
  taus <- sapply(x, function(x) x$tau)
  
  cf = lapply(x, function(x) x$b)
  if (ncol(cf[[1]]) == 4) {
    for (i in 1:length(cf)) {
      cfi <- cf[[i]]
      cfi <- cbind(cfi[, 1], cfi[, 1] - cfi[, 2] * zalpha, 
                   cfi[, 1] + cfi[, 2] * zalpha)
      colnames(cfi) <- c("coefficients", "lower bd", "upper bd")
      cf[[i]] <- cfi
    }
  }
  if (ncol(cf[[1]]) != 3) 
    stop("summary.rqs components have wrong dimension")
  if (is.null(parm)) 
    parm <- rownames(cf[[1]])
  if (is.numeric(parm)) 
    parm <- rownames(cf[[1]])[parm]

  cf <- lapply(cf, function(x) x[parm, , drop = FALSE])
  names(cf) <- paste("tau=", taus)
  
  mfrow_orig <- par("mfrow")
  mar_orig <- par("mar")
  if (is.null(mfrow)) 
    mfrow <- n2mfrow(length(parm))
  if (is.null(mar)) 
    mar <- c(3.1, 3.1, 3.1, 1.6)
  par(mfrow = c(3,3), mar = mar)
  col <- rep(col, length.out = 2)
  lty <- rep(lty, length.out = 2)
  if (is.null(border)) 
    border <- col[2]
  if (is.null(main)) 
    main <- parm

  main <- rep(main, length.out = length(parm))
  xlab <- rep(xlab, length.out = length(parm))
  ylab <- rep(ylab, length.out = length(parm))
  ylim0 <- ylim
  for (i in seq(along = parm)) {
    b <- t(sapply(seq(along = cf), function(tau) cf[[tau]][i, ]))
    if (is.null(ylim))  {
      ylim <- range(b[, 2], b[, 3])
    }
    plot(rep(taus, 2), c(b[, 2], b[, 3]), type = "n", ylim = ylim, 
         xlab = xlab[i], ylab = ylab[i], main = main[i])
    polygon(c(taus, rev(taus)), c(b[, 2], rev(b[, 3])), col = col[2], 
            border = border)
    points(taus, b[, 1], cex = cex, pch = pch, type = type, 
           col = col[1])
    abline(h = 0, col = gray(0.3))
    ylim <- ylim0
  }
  par(mfrow = mfrow_orig, mar = mar_orig)
  x <- cf
  invisible(structure(as.vector(unlist(x)), .Dim = c(dim(x[[1]]), length(x)), .Dimnames = list(rownames(x[[1]]), colnames(x[[1]]), 
                                                                                               names(x))))
}

# Later we may consider the normalization of the factor and factor loading estimation???
PC_est<-function(Res,R,n,t){
  # need to be modified to select the number of factors automatically, which is however not a big issue
  Z<-matrix(0,t,n) # T by N matirx and therefore each colulmn is data for an individual
  for (i in 1:n){
    Z[,i]<-Res[(t*(i-1)+1):(t*i)]
  }
  VEC <- eigen(Z%*%t(Z))$vectors; 
  F <- sqrt(t)*(VEC)[,1:R]; # The factors is T by r matrix
  #L <- t(solve(t(F)%*%F)%*%t(F)%*%Z); # The factor loading is are N by r matrix
  #L <- t(Z)%*%F/t; # follow bai and Ng 2002
  L <- t(solve(t(F)%*%F)%*%t(F)%*%Z);#Follow Ando and Bai (2019)
  FL <- F%*%t(L)
  list(FL=FL,F=F,L=L)
}

Mean_IE<-function(X_final,AY,R,iter,cri_V){
  
    # factors are quantile invariant (i.e not quantile dependent)
    # X,y are NT by k and (T+1) by N panel data with N individuals and T time periods;
    # where (T+1) is due to the fact that lagged term consumes one degree of freedom
    # tau is the quatile; R is the number of factors (later we may consider estimated version)
    # iter is the number of iteration times for IPC; cri_v is the creterion for quitting iteration
    ##First Step: iterated estimation of factors
    #1. Initial value
    n<-ncol(AY)
    t<-nrow(AY)-1 # minus 1 is for lag
    y <- vec(AY[2:nrow(AY),])
    X <-cbind(vec(AY[1:(nrow(AY)-1),]), X_final)
    fit_initial<- lm(y~X)
    B_initial<- fit_initial$coefficients
    res_initial<-fit_initial$residuals
    
    
    FL <- PC_est(res_initial,R,n,t)$FL
    
    B.old <- B_initial
    FL.old <- FL
    
    for(ITE in 1:iter){
      fit <- lm(y-vec(FL)~X) #plus 0 means no constant in the model, but here we need the constant
      B <- fit$coefficients
      #XB<- X%*%B
      XB<- cbind(1,X)%*%B
      res <- y-XB
      FL <- PC_est(res,R,n,t)$FL
      F<-PC_est(res,R,n,t)$F
      L_PC<-PC_est(res,R,n,t)$L # We also need to store the estimation of factor loadings and coefficient 
      B_PC<-B # for bias correction and calculation of variance.
      e_ols<-fit$residuals
      #Up <- sum((B.old-B)^2)/(sum(B.old^2))+sum((FL.old-FL)^2)/(sum(FL.old^2)) 
      #Up <- sum((B.old-B)^2)/(length(B))+sum((FL.old-FL)^2)/(length(FL))
      Up <- mean((B.old-B)^2)+mean((FL.old-FL)^2)
      
      B.old <- B
      FL.old <-FL
      
      #print(Up)
      #print(ITE)
      
      if(Up<=cri_V){break}
    }
    
  ICvalue<-matrix(0,1,3)
  
  ICvalue[,1]=log(mean(e_ols^2))+R*(n+t)*log(n*t/(n+t))/(n*t)
  ICvalue[,2]=log(mean(e_ols^2))+R*(n+t)*log(min(n,t))/(n*t)
  ICvalue[,3]=log(mean(e_ols^2))+R*log(min(n,t))/(min(n,t))
  
  #Store all the estimation results
  list(IC_val=ICvalue,B=B_PC)
}




Mean_IE2<-function(X_final,AY,R,iter,cri_V){
  
  # factors are quantile invariant (i.e not quantile dependent)
  # X,y are NT by k and (T+1) by N panel data with N individuals and T time periods;
  # where (T+1) is due to the fact that lagged term consumes one degree of freedom
  # tau is the quatile; R is the number of factors (later we may consider estimated version)
  # iter is the number of iteration times for IPC; cri_v is the creterion for quitting iteration
  ##First Step: iterated estimation of factors
  #1. Initial value
  n<-ncol(AY)
  t<-nrow(AY)-1 # minus 1 is for lag
  FL <- PC_est(AY[2:nrow(AY),],R,n,t)$FL
  
  y <- vec(AY[2:nrow(AY),])
  X <-cbind(vec(AY[1:(nrow(AY)-1),]), X_final)
  fit_initial<- lm(y-vec(FL)~X)
  B_initial<- fit_initial$coefficients
  res_initial<-fit_initial$residuals
  
  
  B.old <- B_initial
  FL.old <- FL
  
  for(ITE in 1:iter){
    XB<- cbind(1,X)%*%B.old
    res <- y-XB
    FL <- PC_est(res,R,n,t)$FL
    fit <- lm(y-vec(FL)~X) #plus 0 means no constant in the model, but here we need the constant
    B <- fit$coefficients
    #XB<- X%*%B
    # XB<- cbind(1,X)%*%B
    # res <- y-XB
    # FL <- PC_est(res,R,n,t)$FL
    # F<-PC_est(res,R,n,t)$F
    # L_PC<-PC_est(res,R,n,t)$L # We also need to store the estimation of factor loadings and coefficient 
    # B_PC<-B # for bias correction and calculation of variance.
    # e_ols<-fit$residuals
    #Up <- sum((B.old-B)^2)/(sum(B.old^2))+sum((FL.old-FL)^2)/(sum(FL.old^2)) 
    #Up <- sum((B.old-B)^2)/(length(B))+sum((FL.old-FL)^2)/(length(FL))
    Up <- mean((B.old-B)^2)+mean((FL.old-FL)^2)
    
    B.old <- B
    FL.old <-FL
    
    #print(Up)
    #print(ITE)
    
    if(Up<=cri_V){break}
  }
  
  #Store all the estimation results
  list(B=B)
}





rq.fit.PanelIE <- function(X_final,AY,tau=TAU,R,iter=200,cri_V=0.00001,c1,power1,smooth,dynamic){
  # factors are quantile invariant (i.e not quantile dependent)
  # X,y are NT by k and (T+1) by N panel data with N individuals and T time periods;
  # where (T+1) is due to the fact that lagged term consumes one degree of freedom
  # tau is the quatile; R is the number of factors (later we may consider estimated version)
  # iter is the number of iteration times for IPC; cri_v is the creterion for quitting iteration
  ##First Step: iterated estimation of factors
  #1. Initial value
  n<-ncol(AY)
  t<-nrow(AY)-1 # minus 1 is for lag
  y <- vec(AY[2:nrow(AY),])
  if(dynamic=="T"){
    X <-cbind(vec(AY[1:(nrow(AY)-1),]), X_final)
  } else{
    X <-X_final
  }

  fit_initial<- lm(y~X)
  B_initial<- fit_initial$coefficients
  res_initial<-fit_initial$residuals
  
  
  FL <- PC_est(res_initial,R,n,t)$FL
  
  B.old <- B_initial
  FL.old <- FL
  
  for(ITE in 1:iter){
    fit <- lm(y-vec(FL)~X) #plus 0 means no constant in the model, but here we need the constant
    B <- fit$coefficients
    #XB<- X%*%B
    XB<- cbind(1,X)%*%B
    res <- y-XB
    FL <- PC_est(res,R,n,t)$FL
    F<-PC_est(res,R,n,t)$F
    L_PC<-PC_est(res,R,n,t)$L # We also need to store the estimation of factor loadings and coefficient 
    B_PC<-B # for bias correction and calculation of variance.
    e_ols<-fit$residuals
    #Up <- sum((B.old-B)^2)/(sum(B.old^2))+sum((FL.old-FL)^2)/(sum(FL.old^2)) 
    #Up <- sum((B.old-B)^2)/(length(B))+sum((FL.old-FL)^2)/(length(FL))
    Up <- mean((B.old-B)^2)+mean((FL.old-FL)^2)
    
    B.old <- B
    FL.old <-FL
    
    #print(Up)
    #print(ITE)
    
    if(Up<=cri_V){break}
  }
  
  ###Second, substitute the factors and run quantile regression
  
  IN<-diag(1,n,n)
  IT<-diag(1,t,t)
  IF<-kronecker(IN, F)
  
  #Alternative way (This seems to be the best way)
  #Instead of defactor first, we directly write the model using dummy variable
  #We therefore could obtain quantile dependent coefficient for beta
  #quantile and individual specific coefficient for factor loadings
  y_final<- vec(AY[2:nrow(AY),])
  if(dynamic=="T"){
    X_final<-cbind(vec(AY[1:(nrow(AY)-1),]), X_final,IF) 
  } else{
    X_final<-cbind(X_final,IF)
  }
  
  fit_final <- rq(y_final~X_final,TAU)
  
  
  
  
  if(smooth == "Ture"){
    s1<-var(fit_final$residuals)
    ##Following Galvao(2016), the bandwidth dependens on the variance of residual from non-smoothed regression
    h<-c1*s1*(n*t)^(power1)
    # Smoothed estimator
    # Written according to Shi(Zhentao) (2019) and Chen (2019)
    X_smooth<-cbind(1,X_final)
    ## Basically, the following function calculate the loss of the objective function at a specific parameter value
    ## like the likelihood/ sum square of residual etc.
    ##The calculation is based on Chen(2019a,b) and Kato and Galvao(2016) and chen and Huo (2019)
    qr.Smooth<-function(b) {
      b = as.matrix( b )
      z = (y_final- X_smooth %*% b)/h
      # we may also consider the following kernel function  
      #3465/8192 *(7-105*z^2+462*z^4-858*z^6+715*z^8-221*z^10)
      ##Integration function
      k <- function (t) {
        ifelse(abs(t)<=1,1,0)*3465/8192 *(7-105*t^2+462*t^4-858*t^6+715*t^8-221*t^10) #we need eigth order kernel function
        ###105/64 *(1-5*t^2+7*t^4-3*t^6) 
      }
      kz<-matrix(0,n*t,1)
      for (i in 1:(n*t)) {kz[i,]<-1-integrate(k,  lower=-1, upper=z[i,])$value}
      qr.smooth.obj<-(TAU-kz)*z*h ##It is importance to notice the definition of z
      sum.qr.smooth.obj = sum(qr.smooth.obj)
      return(sum.qr.smooth.obj)
    }
    
    #b.init = c(0.1,1,2,2,0.5,0.01) # initial value
    #Following the STARF project, we may rewrite the folloiwng code using alternative methods and parallel optimization
    b.init<-fit_final$coefficients
    b.hat = optimx(b.init, qr.Smooth, method = c("Nelder-Mead"), control = list(reltol = 1e-7, abstol = 1e-7))
    
    list(fit=fit_final,sfit=b.hat,B_OLS=B_PC,L_OLS=L_PC,F=F,e_ols=e_ols)
  }
  else{
    list(fit=fit_final,B_OLS=B_PC,L_OLS=L_PC,F=F,e_ols=e_ols) 
  }
}

#### Function for calculating the norm
norm<-function(x){sqrt(sum(x^2))}


k <- function (t) {
  ifelse(abs(t)<=1,1,0)*3465/8192 *(7-105*t^2+462*t^4-858*t^6+715*t^8-221*t^10) #we need eigth order kernel function
  ###105/64 *(1-5*t^2+7*t^4-3*t^6) 
}

trace_ratio<-function(x,y){
  TR<-tr(t(y)%*%x%*%solve(t(x)%*%x)%*%t(x)%*%y)/tr(t(y)%*%y)
  return(TR)
}


### For calculation of CSA
PQ <- function(h, id){
  if(is.vector(h))
  {h <- matrix(h, ncol = 1)}
  Ph <- unique(id)
  Ph <- cbind(Ph, table(id))
  for(i in 1:ncol(h))
  {Ph <- cbind(Ph, tapply(h[, i], id, mean))}## mean of the same index (cross section average, include constant) 
  #https://stackoverflow.com/questions/3505701/grouping-functions-tapply-by-aggregate-and-the-apply-family
  is <- tapply(id, id)
  Ph <- Ph[is, - (1:2)]
  Qh <- h - Ph ##what is the meaning of Qh
  list(Ph=as.matrix(Ph), Qh=as.matrix(Qh), is=is)##what is the meaning of Qh
}




