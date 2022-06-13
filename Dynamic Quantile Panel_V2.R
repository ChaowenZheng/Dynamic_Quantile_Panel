
####Code for dynamic quantile panel regression 
#### Except factors, all the other parameters are allowed to be quantile dependent
#### Our proposed estimation procedure:
#### Step 1: PC estimation of factors
#### Step 2: substitute these factors back and run quantile regression
## Note: be careful that we allow factor loadings to be both quantile and individual dependent
#### For comparison, we also consider other estimators:
#### 1. Smoothed version of our proposed estimator (Kato and Galvao,2016; chen,2019;Chen and Huo(2019))
#### 2. Bias corrected version using split panel jackknife (Van and Weidner,2016)
#### 3. Standard Quantile Regression that ignores the factors
#### 4. IVQR for dynamic panel with fixed effects (Galvao,2011)
#### 5. PIVQR for dynamic panel with fixed effects (Galvao,2010)
#### 6. Harding, Larmarche and Pesaran (2018), CCE based approach
## (which relies on the assumption that y and x share exactly the same set of factors
## and also regressors are not correlated with factor loadings)
#### 7. Chen (2019), apply PC to Y directly
#### Written by CZ, This version 15/03/2020

rm(list=ls())
# check.packages <- function(pkg){
#   new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#   if (length(new.pkg)) 
#     install.packages(new.pkg, dependencies = TRUE)
#   sapply(pkg, require, character.only = TRUE)
# }

# Usage example
# packages<-c("quantreg", "matrixcalc", "SparseM", "ks", "AER", "numDeriv","optimx","nimble")
# check.packages(packages)

library(quantreg)
#library(quantreg,lib.loc="/path/to/R-packages/")
#https://stackoverflow.com/questions/30034542/unable-to-load-any-package-in-r-unable-to-load-shared-object
library(matrixcalc)#for vectorization
library(SparseM) #for IVQR
library(ks) # for the use of invvec ( inverse operation of vectorization operator)
library(AER) ##The following three package are for nonlinear optimization
library(numDeriv)
library(optimx)
library(psych) # For trace function
library(nimble)

#setwd("//userfs/cz1113/w2k/Desktop/Chapter2")
setwd("Z:/cz1113/Chapter_2/Chapter2_MC/20220129")
#setwd("\\Ecopc1055\c$\Users/cz1113/Google Drive/PhD/Potential Topic/Quantile/Dynamic Quantile Panel Data with Interactive Effects")

#===================================#
##          Define Functions       ##
#===================================#
####First, functions for IVQR and PIVQR
####From Galvao PIVQR.txt
### Problem with this code:
##1. Individual effects is not quantile varying, therefore estimation of several quantiles simulataneously
##2. IVQR is achieved by grid serach for values of autogressive coefficient
##3. How to achieve the above two at the same time is difficult (ongoing challenge)
###The following code achieve 2 but not 1 through estimation only one quantile
###To run Penalized quantile regression, we can either install the following package or use the function give in 
#http://www.econ.uiuc.edu/~roger/research/panel/long.html
#install.packages("rqpd", repos="http://R-Forge.R-project.org")

rq.fit.panel <- function(X,y,s,w,taus,lambda){
  # the matrix X is assumed to contain an intercept
  # the vector s is a strata indicator assumed (so far) to be a one-way layout (individual effects)
  # NB:  
  # 1.  The value of the shrinkage parameter lambda is an open research problem in
  # 	the simplest homogneous settings it should be the ratio of the scale parameters
  # 	of the fixed effects and the idiocyncratic errors
  # 2.  On return the coefficient vector has m*p + n elements where m is the number
  #	quantiles being estimated, p is the number of colums of X, and n is the
  #	number of distinct values of s.  The first m*p coefficients are the 
  #	slope estimates, and the last n are the "fixed effects"
  # 3.  Like all shrinkage (regularization) estimators, asymptotic inference is somewhat
  #	problematic... so the bootstrap is the natural first resort.
  
  require(SparseM)
  require(quantreg)
  K <- length(w)
  if(K != length(taus))
    stop("length of w and taus must match")
  # number of quantiles and number of weights must match
  X <- as.matrix(X)
  p <- ncol(X)
  n <- length(levels(as.factor(s)))
  N <- length(y)
  if(N != length(s) || N != nrow(X))
    stop("dimensions of y,X,s must match")
  Z <- as.matrix.csr(model.matrix(~as.factor(s)-1))
  Fidelity <- cbind(as(w,"matrix.diag.csr") %x% X,w %x% Z)
  Penalty <- cbind(as.matrix.csr(0,n,K*p),lambda*as(n,"matrix.diag.csr"))
  D <- rbind(Fidelity,Penalty)
  y <- c(w %x% y,rep(0,n))
  a <- c((w*(1-taus)) %x% (t(X)%*%rep(1,N)),
         sum(w*(1-taus)) * (t(Z) %*% rep(1,N)) + lambda * rep(1,n))
  rq.fit.sfn(D,y,rhs=a)
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
rq.fit.PanelIE <- function(AX_T,AY,tau=TAU,R,iter=200,cri_V=0.00001,c1,power1,smooth){
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
  X <-cbind(vec(AY[1:(nrow(AY)-1),]), AX_T)
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
  X_final<-cbind(vec(AY[1:(nrow(AY)-1),]), AX_T,IF) 
  
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
####Function for Chen(2019)



rq.fit.PanelIE_Chen <- function(X,AY,tau=TAU,R){
  ###Estimation of Factors
  Cy<-AY[2:nrow(AY),]
  n<-ncol(Cy)
  t<-nrow(Cy)
  CVEC <- eigen(Cy%*%t(Cy))$vectors; 
  CF <- sqrt(T)*(CVEC)[,1:R]; 
  CL <- t(solve(t(CF)%*%CF)%*%t(CF)%*%Cy);
  CFL <- CF%*%t(CL)
  IN<-diag(1,n,n)
  IT<-diag(1,t,t)
  CIF<-kronecker(IN, CF)
  #Estimation of the coefficient
  CX <-cbind(vec(AY[1:(nrow(AY)-1),]),X,CIF)
  fit <- rq(vec(Cy)~CX,TAU)
  list(fit=fit)
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

### Function for analysing whether the factor are well estimated
## See Stock and Watson (2002)
trace_ratio<-function(x,y){
  TR<-tr(t(y)%*%x%*%solve(t(x)%*%x)%*%t(x)%*%y)/tr(t(y)%*%y)
  return(TR)
}

#===================================#
##               Design            ##
#===================================#


#===================================#
##               Design            ##
#===================================#
R <- 2  #number of factors
addr<-2 #number of additional factors
p <- 2  #number of regressors
VAR <- 1  #for variance (may delete later)
a<-0.2 # beta(tau)=beta+a*tau
rho<-0.5 #autoregressive parameter
Xphi<-0.2#autoregressive parameter for error in x 
Fphi<-0.2#autoregressive parameter for factors 
#In Chen (2021), he considers iid factors, here we allow more robust structure.
# we set it as 0.2 currently, but for normal and t, they are Ok even when =0.5
robust<-0
beta<-c(1:p) # slope coefficient for regressors
#T <- 50  #time 
N <- 30 
#individuals
#h<-0.1 #bandwidth of the kernel function used in smoothed quantile estimation
#Following Galvao(2016), we now choose bandwidth according to the data
BN<-20 #length of burn-in time period 
TaB<-T+BN #Total time period
TAU <- 0.2#quantile
MC<-250#Number of MCs
exp<-1 #experiment indicator
dis<-2 #distribution indicator
##1 for normal 2 for t(4), 3 for chi(1), 4 for chi(2) and 5 for chi(3)
## We will not assume a t(2) distribution as the results are not satisfactory.
fl_mean<-0.5
fl_var<-0.5

fl_x_mean<-0.5
fl_x_var<-0.5
fl_cor<-0.8 ##we may consider smaller value

kappa1<-2 #parameter in front of the idiosyncratoc error
kappa2<-2#parameter in front of the idiosyncratoc error
lambda<-0.5 #tuning parameter for PIVQR
TB <- matrix(c(rho,beta)+a*TAU,p+1) #ture parameter matrix
TB_OLS <- matrix(c(rho,beta)+a*0.5,p+1) #ture parameter matrix at mean, where E(u_{it}=0.5)
#This is intuitive and when the error is symetric, the mean un general equals to the 0.5 quantile
#What if the error is not symetrics like chi-square? (think about this more)
#Seems to be also the case since u_{it} is a uniformally distributed variable.
TB<-do.call("cbind", rep(list(TB),MC))
FB <- matrix(c(rho,beta)+a*TAU+0.2,p+1) #ture parameter matrix
FB<-do.call("cbind", rep(list(FB),MC))
iter<-200 # maximum iteration for IPC
cri_V<-10^-7 # Convergence criterion
###Bandwidth selection
###Currently, we do not report Smoothed Results, to save time
###Later we may have a test to see if this is better or not
c1<-2
power1<--1/10


#h_range<-c(3,4,1.5,2,2.5)##this is the bandwidth selection for estimating the sparsity function
#Based on extensive simulations, we finally decide to use on the hf bandwidth with scale as 4
#However, when the error is chi(3), we use scale 2
h_range<-c(3,4)
# we consider runing simulation for multiple time periods to save time
#t_range<-c(rep(30,1))
t_range<-c(30,50,100)
F_hat<-matrix(0,length(t_range),MC)
####define matrix to store the estimation results
SQbeta<-matrix(0,length(t_range)*(p+1),MC)
TSbeta<-matrix(0,length(t_range)*(p+1),MC)
CCEbeta<-matrix(0,length(t_range)*(p+1),MC)
IVQRbeta<-matrix(0,length(t_range)*(p+1),MC)
PIVQRbeta<-matrix(0,length(t_range)*(p+1),MC)
Cbeta<-matrix(0,length(t_range)*(p+1),MC)
beta_c_chen<-matrix(0,length(t_range)*(p+1),MC)
beta_c<-matrix(0,length(t_range)*(p+1),MC)
Smoothbeta<-matrix(0,length(t_range)*(p+1),MC)
##For size performance
t_size_c<-array(0,dim=c(length(t_range)*(p+1),MC,length(h_range)))
t_size_robust<-array(0,dim=c(length(t_range)*(p+1),MC,length(h_range))) 

## For Power performance
t_power_c<-array(0,dim=c(length(t_range)*(p+1),MC,length(h_range)))
t_power_robust<-array(0,dim=c(length(t_range)*(p+1),MC,length(h_range))) 


old_all <- Sys.time() # get start time
t_index<-0
ar<-matrix(0,MC,1)# for storing the estimated ar coefficient using IVQR
for(T in t_range){
  t_index<-t_index+1
  TaB<-T+BN #Total time period
  for(m in 1:MC){
    old <- Sys.time() # get start time
    #===================================#
    ##   Generate simulated Data  ##
    #===================================#
    U <- matrix(runif(TaB*N,0,1),nrow=TaB,ncol=N)  #idiosyncratic error
    
    ##Generate of factors
    
    FAC<-matrix(0,TaB,R+addr)
    for (b in 1:(R+addr)){
      F_e<-matrix(rnorm(TaB,0,(1-Fphi^2)),TaB,1)
      for (t in 2:TaB){
        FAC[t,b]=Fphi*FAC[(t-1),b]+F_e[t,]
      }
    }
    #FAC<-as.matrix(F[(BN+1):nrow(F),])
    
    ###Generate of factor loadings
    LAM<-matrix(rnorm(N*R,fl_mean,fl_var),N,R)# We may assume the mean of factor loadings to be non-zero
    #As a result, the IVQR by Galvao (2010) then could become worse
    # Since the bais correction depends on the estimation of factors at both the mean and quantile level
    # we also need to obtain somewhat consistent estimation for factor loadings for each quantile and individual.
    #While for each individual, we allow the quantiel dependence to be the same for all individuals.
    #We therefore may construct the estimation for factor loadings via \gamma_{i}^{M}+a*tau
    # X1LAM<-matrix(rnorm(N*R,fl_x_mean,fl_x_var),N,R)# pure random factor loadings for x
    # X2LAM<-matrix(rnorm(N*R,fl_x_mean,fl_x_var),N,R) #It is better to assume zero mean for the factor loadings
    X1LAM<-cbind(LAM[,1]*fl_cor,matrix(rnorm(N*1,fl_x_mean,fl_x_var),N,1))
    # allow correlated factor loading in x and y
    X2LAM<-cbind(LAM[,2]*fl_cor,matrix(rnorm(N*1,fl_x_mean,fl_x_var),N,1))
    #It is better to assume zero mean for the factor loadings
    #which give better estimation of factor in the IPC analysis
   
    
    FL <- matrix(0,nrow=TaB,ncol=N) # product of factor and factor loadings
    X1FL <- matrix(0,nrow=TaB,ncol=N) # product of factor and factor loadings for X1
    X2FL <- matrix(0,nrow=TaB,ncol=N) # product of factor and factor loadings for X2
    
    for(t in 1:TaB){
      for(j in 1:N){
        LAM_tau <- LAM[j,]+a*U[t,j] #factor loading is allowed to be quantile varying
        # for each of the factor loading, the quantile dependence is the same
        # but for different time, the data can be seen as generated at different quantiles (u_{it}) 
        FL[t,j] <- FAC[t,1:R]%*%LAM_tau
        
        
        ##for X
        if(R ==2){
          X1FL[t,j]<- FAC[t,c(1,3)]%*%X1LAM[j,]
          X2FL[t,j]<- FAC[t,c(2,4)]%*%X2LAM[j,]
        }
        if(R ==3){
          X1FL[t,j]<- FAC[t,c(1,3,4)]%*%X1LAM[j,]
          X2FL[t,j]<- FAC[t,c(1,2,5)]%*%X2LAM[j,]
        }
      }
    }
    
    ############ Generate independent variable
    ##Step 1: generate error
    #We allow the errors to be serial correlated
    
    AX<-matrix(0,TaB*N,p)
    AX_tem<-matrix(0,TaB,p)
    for (i in 1:N){
      X_e<-matrix(rnorm(TaB*p,0,(1-Xphi^2)),TaB,p)
      for (t in 2:TaB){
        AX_tem[t,]=Xphi*AX_tem[t-1,]+X_e[t,]
      }
      AX[(TaB*(i-1)+1):(TaB*i),]<-AX_tem
    }
    
    ### Step 2: allow X to be correlated with factors and factor loadings
    
    for(i in 1:N){
      for(t in 1:TaB){
        # xf <- FAC[t,1]^2; xl <- LAM[i,1]^2
        # AX[TaB*(i-1)+t,1] <- AX[TaB*(i-1)+t,1]+0.02*xf+0.02*xl
        # xf <- FAC[t,2]^2; xl <- LAM[i,2]^2
        # AX[TaB*(i-1)+t,2] <- AX[TaB*(i-1)+t,2]-0.01*xf+0.02*xl
        #the factor loadings of x on factors play a role. For more discussion, see Westlund (2015)
        # For a better DGP specification, we may also refer to KSS (2020), test for factor loading correlation
        #xf1 <- FAC[t,1]; xl1 <- LAM[i,1]; xf2 <- FAC[t,2]; xl2 <- LAM[i,2]
        # AX[TaB*(i-1)+t,1] <- 0.5*AX[TaB*(i-1)+t,1]+0.2*xf1+1*xl1
        # AX[TaB*(i-1)+t,2] <- 0.5*AX[TaB*(i-1)+t,2]-0.2*xf2+1*xl2
        ###According to Westerlund(2013), it seem that the following way of generating correlated factor loadings is problematic
        ### One the one hand, it weakens the endogeneity issues acused by X
        ### On the other hand, correlated factor loadings means factor loadings are correlated, not means X are correlated with factor loading.
        AX[TaB*(i-1)+t,1] <- kappa2*AX[TaB*(i-1)+t,1]+X1FL[t,i]#+X12LAM[i,]*LAM[i,1]
        AX[TaB*(i-1)+t,2] <- kappa2*AX[TaB*(i-1)+t,2]+X2FL[t,i]#+X22LAM[i,]*LAM[i,2]
      }
    }
    AX_T<-matrix(0,T*N,p)
    for(i in 1:N){
      AX_T[(T*(i-1)+1):(T*i),]<- AX[(TaB*(i-1)+BN+1):(TaB*i),]
    }
    
    ######Generate the product of XB
    XB <- matrix(0,nrow=TaB,ncol=N)
    
    
    for(i in 1:N){
      X <- AX[(TaB*(i-1)+1):(TaB*i),]
      for(t in 1:TaB){
        B <-beta+a*U[t,i]
        #later we may consider heterogeneous coefficient as follows
        #B <- c(1,2)+i/N+a*U[t,i]
        XB[t,i] <- X[t,]%*%B
      }
    }
    
    ###later we may consider heterogeneous case
    #TB <- matrix(0,p,N)
    #for(j in 1:N){TB[,j] <- c(1,2)+j/N+a*TAU}
    
    ###Generate of Y
    
    AY<-matrix(0,TaB,N)
    for (t in 2:TaB){
      for (i in 1:N){
        if(dis==1){
          AY[t,i] <- (rho+a*U[t,i])*AY[t-1,i]+XB[t,i]+FL[t,i]+kappa1*qnorm(U[t,i],0,VAR) ##include constant
        }
        ## t distribution (It is interesting to see that results for t(2) are strange)
        if(dis==2){
          AY[t,i] <- (rho+a*U[t,i])*AY[t-1,i]+XB[t,i]+FL[t,i]+kappa1*qt(U[t,i],4) ##include constant
        }
        ## chi square distribution
        if(dis==3){
          AY[t,i] <- (rho+a*U[t,i])*AY[t-1,i]+XB[t,i]+FL[t,i]+kappa1*(qchisq(U[t,i],1)-1)/sqrt(2) ##include constant demean
        }
        if(dis==4){# This lognormal distribution also does not work.
          AY[t,i] <- (rho+a*U[t,i])*AY[t-1,i]+XB[t,i]+FL[t,i]+kappa1*(qchisq(U[t,i],2)-2)/2 ##include constant
        }
        if(dis==5){# This lognormal distribution also does not work.
          AY[t,i] <- (rho+a*U[t,i])*AY[t-1,i]+XB[t,i]+FL[t,i]+kappa1*(qchisq(U[t,i],3)-3)/sqrt(6) 
        }
        
        
        
        ## chi square distribution
        #AY[t,i] <- (rho+a*U[t,i])*AY[t-1,i]+XB[t,i]+FL[t,i]+qchisq(U[t,i],1) ##include constant
        #Double exponential (acceptable)
        #AY[t,i] <- (rho+a*U[t,i])*AY[t-1,i]+XB[t,i]+FL[t,i]+qdexp(U[t,i], location = 0, scale = 1)##include constant
        # cauchy distribution
        #AY[t,i] <- (rho+a*U[t,i])*AY[t-1,i]+XB[t,i]+FL[t,i]+qcauchy(U[t,i], location = 0, scale = 1)##include constant
        #Beta distribution (does not work)
        #AY[t,i] <- (rho+a*U[t,i])*AY[t-1,i]+XB[t,i]+FL[t,i]+qbeta(U[t,i], 0.5,0.5)-0.5##include constant
        ##lognormal distribution
        #AY[t,i] <- (rho+a*U[t,i])*AY[t-1,i]+XB[t,i]+FL[t,i]+qlnorm(U[t,i], meanlog = 0, sdlog = 1) ##include constant
        
      }
    }

    AY<-AY[BN:nrow(AY),]
    
    #===================================#
    #Estimation 
    #===================================#
    ########Standard QR estimator
    
    y <- vec(AY[2:nrow(AY),])
    X <-cbind(vec(AY[1:(nrow(AY)-1),]),AX_T)
    fit <- rq(y~X,TAU)## constant should be included
    SQb <- fit$coefficients[2:(p+2)]
    SQbeta[((t_index-1)*(p+1)+1):(t_index*(p+1)),m]<-SQb
    
    ########## IVQR and PIVQR
    y <- vec(AY[2:nrow(AY),])
    x <-cbind(1, AX_T)
    y.lag1<-vec(AY[1:(nrow(AY)-1),])
    #Here we use X_{i,t-1} as instruments
    iv <-matrix(0,T*N,p)
    for(i in 1:N){
      iv[(T*(i-1)+1):(T*i),]<- AX[(TaB*(i-1)+BN):((TaB*i)-1),]
    }
    #IV and exo variable
    Y.lag_IV <- cbind(iv,x)
    taus<-TAU
    #Generate of individual effects
    s <- rep(1:N,rep(T,N)) ## rep 1:n and for each, we rep it rep(t,n) times
    
    
    AA<-seq(0.45,0.75,0.005)
    iv_coef<-matrix(0,length(AA),p)
    g<-matrix(0, length(AA),1)
    
    ####IVQR
    for (k in 1:length(AA))
      iv_coef[k,]<-rq.fit.panel(Y.lag_IV,(y-AA[k]*y.lag1),s,1,TAU,0)$coef[c(1:p)]
    end
    
    for (k in 1:length(AA))
      g[k,]<-norm(iv_coef[k])
    end
    
    I<-which.min(g)
    ar[m,]<-AA[I]
    #plot(AA,g)
    
    
    final.reg<-rq.fit.panel(Y.lag_IV,(y-ar[m,]*y.lag1),s,1,TAU,0)#No penalty
    IVQRb<-c(ar[m,],final.reg$coef[4:5])# the first two are coefs for iv; the third is the intercept
    IVQRbeta[((t_index-1)*(p+1)+1):(t_index*(p+1)),m]<-IVQRb
    
    ####PIVQR
    for (k in 1:length(AA))
      iv_coef[k,]<-rq.fit.panel(Y.lag_IV,(y-AA[k]*y.lag1),s,1,TAU,lambda)$coef[c(1:p)]
    end
    
    for (k in 1:length(AA))
      g[k,]<-norm(iv_coef[k])
    end
    
    I<-which.min(g)
    Par<-AA[I] #Penalized ar coef
    #plot(AA,g)
    
    
    Pfinal.reg<-rq.fit.panel(Y.lag_IV,(y-Par*y.lag1),s,1,TAU,lambda)
    PIVQRb<-c(Par,Pfinal.reg$coef[4:5])
    PIVQRbeta[((t_index-1)*(p+1)+1):(t_index*(p+1)),m]<-PIVQRb
    
    ########   Our Proposed estimator
    
    ###  First, run IPC to obtain estimator for factors
    #In Galvao(2016), the bandwidth is choosen based on the variance of residual of non-smoothed quantile regression
    ##For bandwidth in the smoothed regression
    c1<-2
    power1<--1/10
    fit_all<-rq.fit.PanelIE(AX_T,AY,TAU,R,iter,cri_V,c1,power1, smooth="F")
    fit<-fit_all$fit
    #b_smooth<-t(as.matrix(fit_all$sfit[2:4]))
    b <- fit$coefficients[1:4]
    TSbeta[((t_index-1)*(p+1)+1):(t_index*(p+1)),m]<-b[2:4]
    #Smoothbeta[((t_index-1)*(p+1)+1):(t_index*(p+1)),m]<-b_smooth
    
    F<-fit_all$F
    F_hat[t_index,m]<-trace_ratio(F,FAC[(BN+1):nrow(FAC),1:R])
    ### Bias Corrected Estimator Using Split Panel Jackknife
    # The definition can be found in Chen (2019), equation (9)
    
    X1<-invvec(AX_T[,1],nrow = T,ncol=N)
    X2<-invvec(AX_T[,p],nrow = T,ncol=N)
    ###Correction along the cross section dimension
    AY_1N<-AY[,1:(N/2)]
    X1_1N<-X1[,1:(N/2)]
    X2_1N<-X2[,1:(N/2)]
    X_1N<-cbind(vec(X1_1N),vec(X2_1N))
    fit_1N_all<-rq.fit.PanelIE(X_1N,AY_1N,TAU,R,iter,cri_V,c1,power1,smooth="F")
    ##Should smooth be true in this context???
    fit_1N<-fit_1N_all$fit
    B_1N<- fit_1N$coefficients[1:4]
    
    
    
    AY_2N<-AY[,((N/2)+1):N]
    X1_2N<-X1[,((N/2)+1):N]
    X2_2N<-X2[,((N/2)+1):N]
    X_2N<-cbind(vec(X1_2N),vec(X2_2N))
    fit_2N_all<-rq.fit.PanelIE(X_2N,AY_2N,TAU,R,iter,cri_V,c1,power1,smooth="F")
    fit_2N<-fit_2N_all$fit
    B_2N<- fit_2N$coefficients[1:4]
    
    ###Correction along the time dimension
    AY_1T<-AY[1:((T/2)+1),]
    X1_1T<-X1[1:(T/2),]
    X2_1T<-X2[1:(T/2),]
    
    
    X_1T<-cbind(vec(X1_1T),vec(X2_1T))
    fit_1T_all<-rq.fit.PanelIE(X_1T,AY_1T,TAU,R,iter,cri_V,c1,power1,smooth="F")
    fit_1T<-fit_1T_all$fit
    B_1T<- fit_1T$coefficients[1:4]
    
    
    
    AY_2T<-AY[((T/2)+1):(T+1),]
    X1_2T<-X1[((T/2)+1):T,]
    X2_2T<-X2[((T/2)+1):T,]
    X_2T<-cbind(vec(X1_2T),vec(X2_2T))
    fit_2T_all<-rq.fit.PanelIE(X_2T,AY_2T,TAU,R,iter,cri_V,c1,power1,smooth="F")
    fit_2T<-fit_2T_all$fit
    B_2T<- fit_2T$coefficients[1:4]
    
    ## Bias for bias corrected estimator
    
    beta_spjc<-3*b-0.5*(B_1N+B_2N)-0.5*(B_1T+B_2T)
    beta_c[((t_index-1)*(p+1)+1):(t_index*(p+1)),m]<-beta_spjc[2:4]
    
   
    ##########Estimation of Variance
    #The estimation of residual and factor loadings
    residual_tau<-matrix(0,T,N)

    loading_tau<-matrix(fit$coefficients[5:(N*R+4)],N,R,byrow = T)
    residual_tau<-matrix(fit$residuals,T,N)
     
    #Bandwidth estimation for estimating the density function
        h_index<-0
        for(h_scale in h_range){
          h_index<-h_index+1
          if(h_index<=2){
            hf<-T^(-1/3)*(qnorm(0.975))^(2/3)*(1.5*(dnorm(qnorm(TAU))^2)/(2*(qnorm(TAU))^2+1))^(1/3)
            if(dis==5){
              hf<-2*hf
            } else{
              hf<-h_scale*hf
            }
          } else {
            hf<-T^(-1/5)*(4.5*(dnorm(qnorm(TAU))^4)/((2*(qnorm(TAU))^2+1)^2))^(1/5)
            hf<-h_scale*hf
          }
         
          

          ###Estimation of w_{it}
          ###It seems that the estimation of w_{it} is free of sparse density function
          ###But we still follow Galvao(2016) which involves the sparse density function
          k <- function (t) {
            ifelse(abs(t)<=1,1,0)*3465/8192 *(7-105*t^2+462*t^4-858*t^6+715*t^8-221*t^10) #we need eigth order kernel function
            ###105/64 *(1-5*t^2+7*t^4-3*t^6) 
          }
          F<-fit_all$F
          X <-cbind(1,vec(AY[1:(nrow(AY)-1),]),AX_T)
          W<-matrix(0,T*N,p+2)
          for(i in 1:N){
            Xi_i<-matrix(0,p+2,R)
            Omega_i<-matrix(0,R,R)
            for (t in 1:T){
                Xi_i<-Xi_i+k(residual_tau[t,i]/hf)/hf*as.matrix(X[(i-1)*T+t,])%*%F[t,]
                Omega_i<-Omega_i+k(residual_tau[t,i]/hf)/hf*as.matrix(F[t,])%*%F[t,]
                ###Following Galvao, we may only keep k(residual_tau[t,i]/hf)/hf that is larger than 0.01
            }
            for(t in 1:T){
              W[(i-1)*T+t,]<-as.matrix(X[(i-1)*T+t,])-Xi_i%*%solve(Omega_i)%*%F[t,]
            }
          }
          ###Estimation of Delta
          Delta<-matrix(0,p+2,p+2)
          for (i in 1:N) {
              for(t in 1:T){
                Delta<-Delta+k(residual_tau[t,i]/hf)/hf*as.matrix(W[(i-1)*T+t,])%*%W[(i-1)*T+t,]
              }
          }
          Delta<-Delta/(N*T)
          if(is.singular.matrix(Delta) == T){
            diag(Delta)<-diag(Delta)+0.01
          }
          
          
          
          ##########################################################################
          # Estimation of the variance
          # According to Theorem 2, the variance is composed of three terms
          # Moreover, both this two terms are independently distributed across i and t (may not)
          # We therefore consider estimating the varice by treating the two as one term like \tidle{e}=\mu+e
          # The variance is then \sum_{i=1}^N\sum_{t}^T \tidle{e}_{it}^2 since it has zero mean and iid
          
          ###First, compute the \tidle{e}=var_Mtheta+\mu+e
          ###Computation of \mu
          mu_hat<-matrix(0,T*N,p+2)
          for (i in 1:N){
            for (t in 1:T){
              mu_hat[(i-1)*T+t,]<-(TAU-ifelse(residual_tau[t,i]<=0,1,0))*W[(i-1)*T+t,]
            }
          }
          mean(mu_hat)
          
          ###Computation of e
          #OLS factor loading
          L_ols<-fit_all$L_OLS
          e_ols<-fit_all$e_ols                                                                                                                                                                                                                                                                                                                                                                                                    
          e_hat<-matrix(0,T*N,p+2)
          
          for (t in 1:T){
            Lambda_t<-matrix(0,p+2,R)
            for (j in 1:N) {
                Lambda_t<-Lambda_t+k(residual_tau[t,j]/hf)/hf*as.matrix(W[(j-1)*T+t,])%*%loading_tau[j,]
            }
            Lambda_t<-Lambda_t/N
            for (i in 1:N){
              e_hat[(i-1)*T+t,]<-Lambda_t%*%solve(t(L_ols)%*%L_ols/N)%*%as.matrix(L_ols[i,])*e_ols[(i-1)*T+t]
            }
          }
          ##(2021sep)
          ###Computation of var_Mtheta
          #Calculation of D(z)
          
          #MF_OLS
          I_T<-diag(1,T,T)
          MF<-I_T-F%*%solve(t(F)%*%F)%*%t(F)
          D_z<-matrix(0,p+2,p+2)
          for ( i in 1:N){
            tidle_z_i<-MF%*%X[((i-1)*T+1):(i*T),]
            for(j in 1:N){
              Lambda_j_ols<-(t(L_ols[j,]%*%solve(t(L_ols)%*%L_ols/N)%*%L_ols[j,]))[1,1]
              tidle_z_i<-tidle_z_i-Lambda_j_ols*MF%*%X[((j-1)*T+1):(j*T),]/N
            }
            D_z<-D_z+t(tidle_z_i)%*%tidle_z_i
          }
          D_z<-D_z/(N*T)
          #Claculation of \Phi
          
          
          Phi<-matrix(0,p+2,p+2)
          for(t in 1:T){
            Lambda_t<-matrix(0,p+2,R)
            for (i in 1:N) {
                Lambda_t<-Lambda_t+k(residual_tau[t,i]/hf)/hf*as.matrix(W[(i-1)*T+t,])%*%loading_tau[i,]
            }
            Lambda_t<-Lambda_t/N
            for (i in 1:N){
              Phi<-Phi+Lambda_t%*%solve(t(L_ols)%*%L_ols/N)%*%as.matrix(L_ols[i,])%*%X[(i-1)*T+t,]
            }
          }
          Phi<-Phi/(N*T)
          var_Mtheta<-matrix(0,T*N,p+2)
          for (i in 1:N){
            tidle_z_i<-MF%*%X[((i-1)*T+1):(i*T),]
            for(j in 1:N){
              Lambda_j_ols<-(t(L_ols[j,]%*%solve(t(L_ols)%*%L_ols/N)%*%L_ols[j,]))[1,1]
              tidle_z_i<-tidle_z_i-Lambda_j_ols*MF%*%X[((j-1)*T+1):(j*T),]/N
            }
            for(t in 1:T){
              var_Mtheta[(i-1)*T+t,]<-Phi%*%solve(D_z)%*%as.matrix(tidle_z_i[t,])*e_ols[(i-1)*T+t]
            }
          }
          
          
          # ###Computation of sum and variance (two terms)
          # variance_tau<-matrix(0,p+2,p+2)
          # for (i in 1:N){
          #   for (t in 1:T){
          #     variance_tau=variance_tau+(e_hat[(i-1)*T+t,]-mu_hat[(i-1)*T+t,])%*%t(e_hat[(i-1)*T+t,]-mu_hat[(i-1)*T+t,])/(N*T)
          #   }
          #   if(robust==1){
          #     #consider serial correlation
          #     for (t in 1:(T-1)){
          #       variance_tau=variance_tau+2*(mu_hat[(i-1)*T+t,]-e_hat[(i-1)*T+t,])%*%t(mu_hat[(i-1)*T+t+1,]-e_hat[(i-1)*T+t+1,])/(N*(T-1))
          #     }
          #   }
          # }
          # variance_tau_2<-solve(Delta)%*%variance_tau%*%solve(Delta)#/(N*T)
          # diag(variance_tau_2)
          # 
          
          ###Computation of sum and variance (three terms)
          variance_tau<-matrix(0,p+2,p+2)
          variance_tau_robust<-matrix(0,p+2,p+2)
          for (i in 1:N){
            for (t in 1:T){
              variance_tau=variance_tau+(-var_Mtheta[(i-1)*T+t,]+e_hat[(i-1)*T+t,]-mu_hat[(i-1)*T+t,])%*%t(-var_Mtheta[(i-1)*T+t,]+e_hat[(i-1)*T+t,]-mu_hat[(i-1)*T+t,])/(N*T)
              variance_tau_robust=variance_tau_robust+(-var_Mtheta[(i-1)*T+t,]+e_hat[(i-1)*T+t,]-mu_hat[(i-1)*T+t,])%*%t(-var_Mtheta[(i-1)*T+t,]+e_hat[(i-1)*T+t,]-mu_hat[(i-1)*T+t,])/(N*T)
            }
            #if(robust==1){
              #consider serial correlation
              for (t in 1:(T-1)){
                variance_tau_robust=variance_tau_robust+2*(-var_Mtheta[(i-1)*T+t,]+mu_hat[(i-1)*T+t,]-e_hat[(i-1)*T+t,])%*%t(-var_Mtheta[(i-1)*T+t+1,]+mu_hat[(i-1)*T+t+1,]-e_hat[(i-1)*T+t+1,])/(N*(T-1))
              }
            #}
          }
          variance_tau<-solve(Delta)%*%variance_tau%*%solve(Delta)
          variance_tau_robust<-solve(Delta)%*%variance_tau_robust%*%solve(Delta)
        
      
          t_size_c[((t_index-1)*(p+1)+1):(t_index*(p+1)),m,h_index]<-(beta_spjc[2:4]-TB[,m])/sqrt(diag(variance_tau)[2:4])*sqrt(N*T)
          t_size_robust[((t_index-1)*(p+1)+1):(t_index*(p+1)),m,h_index]<-(beta_spjc[2:4]-TB[,m])/sqrt(diag(variance_tau_robust)[2:4])*sqrt(N*T)
          t_power_c[((t_index-1)*(p+1)+1):(t_index*(p+1)),m,h_index]<-(beta_spjc[2:4]-FB[,m])/sqrt(diag(variance_tau)[2:4])*sqrt(N*T)
          t_power_robust[((t_index-1)*(p+1)+1):(t_index*(p+1)),m,h_index]<-(beta_spjc[2:4]-FB[,m])/sqrt(diag(variance_tau_robust)[2:4])*sqrt(N*T)
      }
    
    ##########################################################################
    #######Chen (2019)
    ###  First, run IPC to obtain estimator for factors
    fit_chen<-rq.fit.PanelIE_Chen(AX_T,AY,tau=TAU,R)
    CB <- fit_chen$fit$coefficients[2:4]
    #TAU
    Cbeta[((t_index-1)*(p+1)+1):(t_index*(p+1)),m]<-CB
    
    ###Spli panel Jacknife
    #### Bias Corrected Estimator Using Split Panel Jackknife
    ## The definition can be found in Chen (2019), equation (9)
    
    # dim(AY_1N)
    # dim(A2X<-X_1N)
    # AY<-AY_1N
    fit_1N_chen<-rq.fit.PanelIE_Chen(X_1N,AY_1N,TAU,R)
    fit_1N_chen<-fit_1N_chen$fit
    B_1N_chen<- fit_1N_chen$coefficients[2:4]
    
    fit_2N_chen<-rq.fit.PanelIE_Chen(X_2N,AY_2N,TAU,R)
    fit_2N_chen<-fit_2N_chen$fit
    B_2N_chen<- fit_2N_chen$coefficients[2:4]
    
    
    
    fit_1T_chen<-rq.fit.PanelIE_Chen(X_1T,AY_1T,TAU,R)
    fit_1T_chen<-fit_1T_chen$fit
    B_1T_chen<- fit_1T_chen$coefficients[2:4]
    
    fit_2T_chen<-rq.fit.PanelIE_Chen(X_2T,AY_2T,TAU,R)
    fit_2T_chen<-fit_2T_chen$fit
    B_2T_chen<- fit_2T_chen$coefficients[2:4]
    
    
    ## Bias for bias corrected estimator
    beta_c_chen[((t_index-1)*(p+1)+1):(t_index*(p+1)),m]<-3*CB-0.5*(B_1N_chen+B_2N_chen)-0.5*(B_1T_chen+B_2T_chen)
    
    
    
    ###Harding and Lamarche
    #use cross section average to approximate for factors
    #In Harding, Lamarche and Pesaran (2018), they use (bary, bary_{t-1}, barx)
    
    ###  First, run IPC to obtain estimator for factors
    
    y <- vec(AY[2:nrow(AY),])
    X <-cbind(vec(AY[1:(nrow(AY)-1),]), AX_T)
    
    st <- rep(1:T,N)
    
    F <- PQ(cbind(y,X),st)$Ph##factor approximation
    dim(F)
    dim(cbind(y,X))
    F<-F[1:T,]
    ###print(trace_ratio(F,TFAC))
    
    ###Second, substitute the factors and run quantile regression
    IN<-diag(1,N,N)
    IF<-kronecker(IN, F)
    #######Fit quantile regression
    
    
    y <- vec(AY[2:nrow(AY),])
    X<-cbind(vec(AY[1:(nrow(AY)-1),]), AX_T,IF)
    fit <- rq(y~X,TAU)
    CCEb <- fit$coefficients[2:4]
    CCEbeta[((t_index-1)*(p+1)+1):(t_index*(p+1)),m]<-CCEb
    
    
    cat(sprintf("\"T=%d\" \"m=%d\" \"h=%.1f\" \"tau=%.1f\" \"dis=%.1f\" \n", t_range[t_index],m,h_range[h_index],TAU,dis))
    
    # print elapsed time
    new <- Sys.time() - old # calculate difference
    print(new) # print in nice format
  }
}


TB<-do.call("rbind", rep(list(TB),length(t_range)))
# TB<-rbind(TB,TB,TB)
SQBias<-as.matrix(rowMeans(SQbeta-TB))
Bias<-as.matrix(rowMeans(TSbeta-TB))
SmoothBias<-as.matrix(rowMeans(Smoothbeta-TB))
Bias_c<-as.matrix(rowMeans(beta_c-TB))
CCEBias<-as.matrix(rowMeans(CCEbeta-TB))
IVQRBias<-as.matrix(rowMeans(IVQRbeta-TB))
PIVQRBias<-as.matrix(rowMeans(PIVQRbeta-TB))
Bias_chen<-as.matrix(rowMeans(Cbeta-TB))
Bias_c_chen<-as.matrix(rowMeans(beta_c_chen-TB))

BIAS<-cbind(SQBias, Bias, SmoothBias,Bias_c,CCEBias, IVQRBias,PIVQRBias,Bias_chen,Bias_c_chen)
colnames(BIAS)=c("QR", "Twostep","Smoothed","SPJ","CCE","IVQR","PIVQR","Chen","SPJ_chen")
row.names(BIAS)<-c(rep("T=30",(p+1)),rep("T=50",(p+1)),rep("T=100",(p+1)))
# row.names(BIAS)<-c(rep("T=30",length(t_range)*(p+1)))
# row.names(BIAS)<-c(rep("T=50",length(t_range)*(p+1)))
#row.names(BIAS)<-c(rep("T=100",length(t_range)*(p+1)))

SQRMSE<-sqrt(rowSums((SQbeta-TB)^2)/MC)
RMSE<-sqrt(rowSums((TSbeta-TB)^2)/MC)
SmoothRMSE<-sqrt(rowSums((Smoothbeta-TB)^2)/MC)
RMSE_c<-sqrt(rowSums((beta_c-TB)^2)/MC)

CCERMSE<-sqrt(rowSums((CCEbeta-TB)^2)/MC)
IVQRRMSE<-sqrt(rowSums((IVQRbeta-TB)^2)/MC)
PIVQRRMSE<-sqrt(rowSums((PIVQRbeta-TB)^2)/MC)
RMSE_chen<-sqrt(rowSums((Cbeta-TB)^2)/MC)
RMSE_c_chen<-sqrt(rowSums((beta_c_chen-TB)^2)/MC)

TRMSE<-cbind(SQRMSE, RMSE, SmoothRMSE,RMSE_c,CCERMSE, IVQRRMSE,PIVQRRMSE,RMSE_chen,RMSE_c_chen)
colnames(TRMSE)=c("QR", "Twostep","Smoothed","SPJ","CCE","IVQR","PIVQR","Chen","SPJ_chen")
row.names(TRMSE)<-c(rep("T=30",(p+1)),rep("T=50",(p+1)),rep("T=100",(p+1)))
# row.names(TRMSE)<-c(rep("T=30",length(t_range)*(p+1)))
# row.names(TRMSE)<-c(rep("T=50",length(t_range)*(p+1)))
# row.names(TRMSE)<-c(rep("T=100",length(t_range)*(p+1)))
BIAS
TRMSE
#For size
size<-matrix(0,length(t_range)*(p+1),length(h_range))
size_c<-matrix(0,length(t_range)*(p+1),length(h_range))
size_robust<-matrix(0,length(t_range)*(p+1),length(h_range))
#For power
power<-matrix(0,length(t_range)*(p+1),length(h_range))
power_c<-matrix(0,length(t_range)*(p+1),length(h_range))
power_robust<-matrix(0,length(t_range)*(p+1),length(h_range))
for(h_i in 1:length(h_range)){
  for (i in 1:(length(t_range)*(p+1))) {
    for(j in 1:MC){
      ##for size
      #size[i,h_i]<-ifelse(abs(t_size[i,j,h_i])>=1.96,1,0)+size[i,h_i]
      size_c[i,h_i]<-ifelse(abs(t_size_c[i,j,h_i])>=1.96,1,0)+size_c[i,h_i]
      size_robust[i,h_i]<-ifelse(abs(t_size_robust[i,j,h_i])>=1.95996,1,0)+size_robust[i,h_i]
      ##for power
      #power[i,h_i]<-ifelse(abs(t_power[i,j,h_i])>=1.96,1,0)+power[i,h_i]
      power_c[i,h_i]<-ifelse(abs(t_power_c[i,j,h_i])>=1.96,1,0)+power_c[i,h_i]
      power_robust[i,h_i]<-ifelse(abs(t_power_robust[i,j,h_i])>=1.95996,1,0)+power_robust[i,h_i]
    }
  }
}

# size<-cbind(size/MC,size_c/MC,size_ac/MC)
# row.names(size)<-c(rep("T=30",(p+1)),rep("T=50",(p+1)),rep("T=100",(p+1)))
# colnames(size)<-c(rep("Two-Step",length(h_range)),rep("SPJ",length(h_range)),rep("A&SPJ",length(h_range)))
# size
# power<-cbind(power/MC,power_c/MC,power_ac/MC)
# row.names(power)<-c(rep("T=30",(p+1)),rep("T=50",(p+1)),rep("T=100",(p+1)))
# colnames(power)<-c(rep("Two-Step",length(h_range)),rep("SPJ",length(h_range)),rep("A&SPJ",length(h_range)))
# power

size<-cbind(size_c/MC,size_robust/MC)
#row.names(size)<-c(rep("T=30",(p+1)),rep("T=50",(p+1)),rep("T=100",(p+1)))
colnames(size)<-c(rep("SPJ_c",length(h_range)),rep("SPJ_robust",length(h_range)))
size
power<-cbind(power_c/MC,power_robust/MC)
#row.names(power)<-c(rep("T=30",(p+1)),rep("T=50",(p+1)),rep("T=100",(p+1)))
colnames(power)<-c(rep("SPJ_c",length(h_range)),rep("SPJ_robust",length(h_range)))
power



Bias_RMSE<-rbind(BIAS, TRMSE)
#All_res
Size_Power<-rbind(size, power)

F_Rsquare<-rowMeans(F_hat)
write.csv(F_Rsquare, file= sprintf("Factor,N=%d, T=(%d,%d,%d), R=%d, MC=%d,autocorr=(%.1f,%.1f,%.1f),FL_XY(%.1f,%.1f,%.1f,%.1f,%.1f),kappa=(%.1f,%.1f), (%d).csv",N,t_range[1],t_range[2],t_range[3],R, MC,rho,Fphi,Xphi,fl_mean,fl_x_mean,fl_var,fl_x_var,fl_cor,kappa1,kappa2,exp))
write.csv(Size_Power, file= sprintf("Size & Power,N=%d, T=(%d,%d,%d),h=(%.1f,%.1f,%.1f,%.1f,%.1f,%.1f), R=%d, MC=%d,TAU=%.2f,autocorr=(%.1f,%.1f,%.1f),robust=%d, Distribution=%d,FL_XY(%.1f,%.1f,%.1f,%.1f,%.1f),kappa=(%.1f,%.1f), (%d).csv",N,t_range[1],t_range[2],t_range[3],h_range[1],h_range[2],h_range[3], h_range[4],h_range[5], h_range[6],R, MC,TAU,rho,Fphi,Xphi,robust,dis,fl_mean,fl_x_mean,fl_var,fl_x_var,fl_cor,kappa1,kappa2,exp))
write.csv(Bias_RMSE, file= sprintf("Bias & RMSE, N=%d, T=(%d,%d,%d),h=(%.1f,%.1f,%.1f,%.1f,%.1f,%.1f), R=%d, MC=%d,TAU=%.2f,autocorr=(%.1f,%.1f,%.1f),robust=%d, Distribution=%d,FL_XY(%.1f,%.1f,%.1f,%.1f,%.1f),kappa=(%.1f,%.1f), (%d).csv",N,t_range[1],t_range[2],t_range[3],h_range[1],h_range[2],h_range[3], h_range[4],h_range[5], h_range[6],R, MC,TAU,rho,Fphi,Xphi,robust,dis,fl_mean,fl_x_mean,fl_var,fl_x_var,fl_cor,kappa1,kappa2,exp))

# print elapsed time
new_all <- Sys.time() - old_all # calculate difference
print(new_all) # print in nice format

