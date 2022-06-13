
####Code for dynamic quantile panel regression empircial analysis


rm(list=ls())


#library(quantreg,lib.loc = "H:/rlibs/4.0.2")
.libPaths()
library(quantreg)
#install.packages("quantreg",lib = "C:/Program Files/R/R-4.0.2/library" )
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
library(plm,lib = "c:/apps/r-3.6.1/library" )
library(rlist,lib = "\\\\userfs/cz1113/w2k/R/win-library/3.6")
#setwd("//userfs/cz1113/w2k/Desktop/Chapter2")
#setwd("\\Ecopc1055\c$\Users/cz1113/Google Drive/PhD/Potential Topic/Quantile/Dynamic Quantile Panel Data with Interactive Effects")


#folder = "//Ecopc1055/c$/Users/cz1113/Google Drive/PhD/chapter2_empirical"
#folder = "C:/Users/86188/Google ‘∆∂À”≤≈Ã/PhD/chapter2_empirical"
folder = "Z:/cz1113/chapter2_empirical"

setwd(folder) 
# setwd("//Ecopc1055/c$/Users/cz1113/Google Drive/PhD/chapter2_empirical")
#setwd("C:/Users/86188/Google ‘∆∂À”≤≈Ã/PhD/chapter2_empirical")
#data<-read.csv("eu2020.csv",header = T)
data<-read.csv("oecd2020.csv",header = T)
source("dynamic_quantile_functions.R")
#data<-data[which(data$codeh!="BELLUX"&data$codef!="BELLUX"&data$codeh!="DEU"&data$codef!="DEU"), ]
dim(data)
names(data)
# Number of individuals
N<-length(unique(data$pair_n))
#IN<-diag(1,N,N)
# data$gdppc<-log(data$gdppcf+data$gdppch)
# data$pop<-log(data$popf+data$poph)
# write.csv(data,"eu2020.csv")
#Time periods
T<-max(data$year)-min(data$year)
IT<-diag(1,T,T)
#Trend
trend<-rep(1:(T+1),each=N)
correct<-1
dynamic<-"T"

#Specifying y
#Y<-data$trade
for (Yid in c(3)){
  if(Yid==1){
    Y<-log(data$rexport)
    YD<-"Export"
  }
  if(Yid==2){
    Y<-log(data$rimport)
    YD<-"Import"
  }
  if(Yid==3){
    Y1<-log(data$rexport)
    Y2<-log(data$rimport)
    YD<-"Export & Import"
  }
  if(Yid==4){
    Y<-0.5*(log(data$rexport)+log(data$rimport))
    YD<-"Export & Import_Sum"
  }
  if(Yid==5){
    Y<-(log(data$rexport+data$rimport))
    YD<-"Export & Import_logSum"
  }

  for (Xid in c(1)){
    if (Xid ==1)
    {
      X<-cbind(0.5*(log(data$gdpf)+log(data$gdph)),0.5*(log(data$popf)+log(data$poph)),data$sim,data$distance,data$border,data$language,data$europ,data$ceep)
      XD<-"Sumlog_GDP_Pop & Sim & Dis_Important & Bor & lan & Europ & cee"
      ET<-"True"
    }
    if (Xid ==2)
    {
      #X<-cbind(0.5*(log(data$gdpf)+log(data$gdph)),0.5*(log(data$popf)+log(data$poph)),data$sim,data$distance,data$border,data$language,data$colony,data$europ,data$ceep)
      X<-cbind(0.5*(log(data$gdpf)+log(data$gdph)),data$sim,data$distance,data$border,data$language,data$colony,data$europ,data$ceep)
      XD<-"Sumlog_GD & Sim & Dis_Important & Bor & lan & col & Europ & cee"
      ET<-"True"
    }
    if (Xid ==3)
    {
      X<-cbind(0.5*(log(data$gdpf)+log(data$gdph)),0.5*(log(data$popf)+log(data$poph)),data$sim,data$distance,data$border,data$colony,data$europ,data$ceep)
      XD<-"Sumlog_GDP_Pop & Sim & Dis_Important & Bor & col & Europ & cee"
      ET<-"True"
    }
    
    if (Xid ==4)
    {
      X<-cbind(0.5*(log(data$gdpf)+log(data$gdph)),0.5*(log(data$popf)+log(data$poph)),data$sim,data$rer,data$distance,data$border,data$colony,data$ceep)
      XD<-"Sumlog_GDP_Pop & Sim & rer & Dis_Important & Bor & col & cee"
      ET<-"False"
    }
    
    if (Xid ==5)
    {
      X<-cbind((log(data$gdpf+data$gdph)),(log(data$popf+data$poph)),data$distance,data$border,data$comlang_ethno,data$europ,data$ceep)
      XD<-"lan_ethno & Europ & Dis_Important& logSum_GDP_pop"
      ET<-"True"
    }
    if (Xid ==6)
    {
      X<-cbind((log(data$gdpf+data$gdph)),(log(data$popf+data$poph)),data$sim,data$distance,data$border,data$language,data$europ,data$ceep)
      XD<-"lan & Europ & Dis_Important& logSum_GDP_pop+sim"
      ET<-"True"
    }
    if (Xid ==7)
    {
      X<-cbind((log(data$gdpf+data$gdph)),(log(data$popf+data$poph)),data$sim,data$distance,data$border,data$europ,data$ceep)
      XD<-"No lan & Europ & Dis_Important& logSum_GDP_pop+sim"
      ET<-"True"
    }
    if (Xid ==8)
    {
      X<-cbind((log(data$gdpf+data$gdph)),(log(data$popf+data$poph)),data$sim,data$distance,data$border,data$language,data$colony,data$europ,data$ceep)
      XD<-"lan & Europ & Dis_Important& logSum_GDP_pop+sim+colony"
      ET<-"True"
    }
    if (Xid ==9)
    {
      X<-cbind((log(data$gdpf+data$gdph)),(log(data$popf+data$poph)),data$sim,data$distance,data$border,data$colony,data$europ,data$ceep)
      XD<-"No lan & Europ & Dis_Important& logSum_GDP_pop+sim+colony"
      ET<-"True"
    }
    if (Xid ==10)
    {
      X<-cbind((log(data$gdpf+data$gdph)),(log(data$popf+data$poph)),data$sim,data$rer,data$rfl,data$distance,data$border,data$language,data$europ,data$ceep)
      XD<-"lan & Europ & Dis_Important& logSum_GDP_pop+sim+rer+rfl"
      ET<-"True"
    }
    
    if (Xid ==11)
    {
      X<-cbind((log(data$gdpf+data$gdph)),(log(data$popf+data$poph)),data$sim,data$rer,data$rfl,data$distance,data$border,data$europ,data$ceep)
      XD<-"No lan & Europ & Dis_Important& logSum_GDP_pop+sim+rer+rfl"
      ET<-"True"
    }
    if (Xid ==12)
    {
      X<-cbind((log(data$gdpf+data$gdph)),(log(data$popf+data$poph)),data$sim,data$rfl,data$distance,data$border,data$language,data$europ,data$ceep)
      XD<-"No lan & Europ & Dis_Important& logSum_GDP_pop+sim+rfl"
      ET<-"True"
    }
    
    dim(X)
    p<-ncol(X)
    if(dynamic=="T"){
      pall<-p+2
      padd<-p+3
    } else{
      pall<-p+1
      padd<-p+2
    }
    matplot(X, type = "l")
    
    #fit<-lm(Y~trend)
    if (Yid !=3){
      Y_wide<-t(matrix(Y,N,T+1,))
      Y_wide<-Y_wide[1:(T+1),]
    } else {
      Y1_wide<-t(matrix(Y1,N,T+1,))
      Y1_wide<-Y1_wide[1:(T+1),]
      Y2_wide<-t(matrix(Y2,N,T+1,))
      Y2_wide<-Y2_wide[1:(T+1),]
    }
    
 
    #for level data
    X_wide<-array(0,dim=c(T+1,N,p))
    for (i in 1:p){
      x_tem<-t(matrix(X[,i],N,T+1))
      X_wide[,,i]<-x_tem
    }
    dim(X_wide)
   
    X_wide<-X_wide[2:(T+1),,]
 
    
    X_final<-matrix(0,N*(T),p)
    for (i in 1:p){
      x_tem<-matrix(X_wide[,,i],N*(T),1)
      X_final[,i]<-x_tem
    }
 
    
    iter<-200 # maximum iteration for IPC
    cri_V<-10^-7 # Convergence criterion
    #TAU<-0.5
    
    
    if (Yid ==3){
      Y_wide<-cbind(Y1_wide,Y2_wide)
      X_wide_2<-array(0,dim=c(T,2*N,p))
      for (i in 1:p){
        X_wide_2[,,i]<-cbind(X_wide[,,i],X_wide[,,i])
      }
      X_final<-rbind(X_final,X_final)
      N<-ncol(Y_wide)
      X_wide<-X_wide_2
    }
    
    ##Panel unit root test
# 
#     purtest(Y_wide, pmax = 4, exo = "intercept", test = "madwu")
#     purtest(Y_wide, pmax = 4, exo = "intercept", test = "levinlin")
#     purtest(Y_wide, pmax = 4, exo = "intercept", test = "ips")
#     purtest(Y_wide, pmax = 4, exo = "intercept", test = "Pm")
#     purtest(Y_wide, pmax = 4, exo = "intercept", test = "hadri")
    
    rmax<-6
    ## Determining the number of factors
    IC_value<-matrix(0,rmax,3)
    for (r_id in 1:rmax){
      mean_fit<-Mean_IE(X_final,Y_wide,r_id,iter,cri_V)
      IC_value[r_id,]<-mean_fit$IC_val
      #print(r_id)
    }
    
    r_hat<-matrix(0,1,3)
    for (ic in 1:3){
      r_hat[,ic]<-which.min(IC_value[,ic])
    }
    r_hat
    R<-r_hat[,2]
    #We do not want 1 factor which would lead modification of the code when calculating the variance
    if (R==1){
      R<-2
    }
    # dim(Y_wide)
    # dim(X_final)
    #===================================#
    #Estimation 
    #===================================#
    ########   Our Proposed estimator
    
    ###  First, run IPC to obtain estimator for factors
    #In Galvao(2016), the bandwidth is choosen based on the variance of residual of non-smoothed quantile regression
    ##For bandwidth in the smoothed regression
    taus = seq(0.11, 0.91, 0.1)
    #taus=c(seq(0.11, 0.41, 0.1),0.53,seq(0.61, 0.91, 0.1))
    taus=c(0.5,0.8)
    #taus=c(0.2)
    Netallout = list()
    
    for (q in 1:length(taus)){
      TAU = taus[q]
      cat("Ydata", Yid,  "Xdata", Xid,"Dynamic",dynamic, "Tau", TAU, "\n")   
      c1<-2
      power1<--1/10
      fit_all<-rq.fit.PanelIE(X_final,Y_wide,TAU,R,iter,cri_V,c1,power1, smooth="F",dynamic=dynamic)
      fit<-fit_all$fit
      #b_smooth<-t(as.matrix(fit_all$sfit[2:4]))

        b <- fit$coefficients[1:(pall)]
      b
      #Smoothbeta[((t_index-1)*(p+1)+1):(t_index*(p+1)),m]<-b_smooth
      
      F<-fit_all$F
      
      beta_spjc<-b
      ### Bias Corrected Estimator Using Split Panel Jackknife
      # The definition can be found in Chen (2019), equation (9)
      
      if(correct==1){
      ###Correction along the cross section dimension
      Y_wide_1N<-Y_wide[,1:(N/2)]
      X1_1N<-X_wide[,1:(N/2),1]
      X_1N<-vec(X1_1N)
      for (i in 2:p){
        X_1N_tem<-X_wide[,1:(N/2),i]
        X_1N<-cbind(X_1N,vec(X_1N_tem))
      }
      
      fit_1N_all<-rq.fit.PanelIE(X_1N, Y_wide_1N,TAU,R,iter,cri_V,c1,power1,smooth="F",dynamic=dynamic)
      ##Should smooth be true in this context???
      fit_1N<-fit_1N_all$fit

      B_1N<- fit_1N$coefficients[1:(pall)]
     
      
      
      Y_wide_2N<-Y_wide[,((N/2)+1):N]
      X1_2N<-X_wide[,((N/2)+1):N,1]
      X_2N<-vec(X1_2N)
      for (i in 2:p){
        X_2N_tem<-X_wide[,((N/2)+1):N,i]
        X_2N<-cbind(X_2N,vec(X_2N_tem))
      }
      
      fit_2N_all<-rq.fit.PanelIE(X_2N,Y_wide_2N,TAU,R,iter,cri_V,c1,power1,smooth="F",dynamic=dynamic)
      fit_2N<-fit_2N_all$fit
   
        B_2N<- fit_2N$coefficients[1:(pall)]
     
      
    
      if(ET=="True")
      {
        beta_spjc<-2*b-0.5*(B_1N+B_2N)

        beta_spjc
       } else{
      #   if(Yid==3){
      #     ###Correction along the time dimension
      #     Y_wide_1T<-Y_wide[1:((T/2)+1),]
      #     X1_1T<-cbind(X_wide[1:(T/2),,1],X_wide[1:(T/2),,1])
      #     X_1T<-vec(X1_1T)
      # 
      #     for (i in 2:p){
      #       X_1T_tem<-cbind(X_wide[1:(T/2),,i],X_wide[1:(T/2),,i])
      #       X_1T<-cbind(X_1T,vec(X_1T_tem))
      #     }
      # 
      #     fit_1T_all<-rq.fit.PanelIE(X_1T,Y_wide_1T,TAU,R,iter,cri_V,c1,power1,smooth="F")
      #     fit_1T<-fit_1T_all$fit
      #     B_1T<- fit_1T$coefficients[1:(pall)]
      # 
      # 
      # 
      #     Y_wide_2T<-Y_wide[((T/2)+1):(T+1),]
      #     X1_2T<-cbind(X_wide[((T/2)+1):T,,1],X_wide[((T/2)+1):T,,1])
      #     X_2T<-vec(X1_2T)
      # 
      #     for (i in 2:p){
      #       X_2T_tem<-cbind(X_wide[((T/2)+1):T,,1],X_wide[((T/2)+1):T,,1])
      #       X_2T<-cbind(X_1T,vec(X_2T_tem))
      #     }
      # 
      #     fit_2T_all<-rq.fit.PanelIE(X_2T,Y_wide_2T,TAU,R,iter,cri_V,c1,power1,smooth="F")
      #     fit_2T<-fit_2T_all$fit
      #     B_2T<- fit_2T$coefficients[1:(pall)]
      #   } else{
          ###Correction along the time dimension
          Y_wide_1T<-Y_wide[1:((T/2)+1),]
          X1_1T<-X_wide[1:(T/2),,1]
          X_1T<-vec(X1_1T)
          
          for (i in 2:p){
            X_1T_tem<-X_wide[1:(T/2),,i]
            X_1T<-cbind(X_1T,vec(X_1T_tem))
          }
          
          fit_1T_all<-rq.fit.PanelIE(X_1T,Y_wide_1T,TAU,R,iter,cri_V,c1,power1,smooth="F")
          fit_1T<-fit_1T_all$fit
          B_1T<- fit_1T$coefficients[1:(pall)]
          
          
          
          Y_wide_2T<-Y_wide[((T/2)+1):(T+1),]
          X1_2T<-X_wide[((T/2)+1):T,,1]
          X_2T<-vec(X1_2T)
          
          for (i in 2:p){
            X_2T_tem<-X_wide[((T/2)+1):T,,1]
            X_2T<-cbind(X_1T,vec(X_2T_tem))
          }
          
          fit_2T_all<-rq.fit.PanelIE(X_2T,Y_wide_2T,TAU,R,iter,cri_V,c1,power1,smooth="F")
          fit_2T<-fit_2T_all$fit
          B_2T<- fit_2T$coefficients[1:(pall)]
        #}
        
        
        beta_spjc<-3*b-0.5*(B_1N+B_2N)-0.5*(B_1T+B_2T)
        
        beta_spjc
       }
      }
      
      ## Bias for bias corrected estimator
      
      
      cat("Estimator", b, "\n", "SPJ_Estimator", beta_spjc, "\n")
      
      dis<-2
      ##########Estimation of Variance
      #The estimation of residual and factor loadings
      residual_tau<-matrix(0,T,N)
    
      loading_tau<-matrix(fit$coefficients[(padd):(N*R+pall)],N,R,byrow = T)
      residual_tau<-matrix(fit$residuals,T,N)
      
      ##For size performance
      h_range<-c(3,4)
      t_size_c<-matrix(0,pall,length(h_range))
      t_size_robust<-matrix(0,pall,length(h_range)) 
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
          ifelse(abs(t)<=1,1,0)*3465/8192 *(7-105*t^2+462*t^4-858*t^6+715*t^8-221*t^10) #we need eighth order kernel function
        }
        F<-fit_all$F
        if(dynamic=="T"){
          X <-cbind(1,vec(Y_wide[1:T,]),X_final)
        } else{
          X <-cbind(1,X_final)
        }

        W<-matrix(0,T*N, pall)
        for(i in 1:N){
          Xi_i<-matrix(0,pall,R)
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
        Delta<-matrix(0,pall,pall)
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
        mu_hat<-matrix(0,T*N,pall)
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
        e_hat<-matrix(0,T*N,pall)
        
        for (t in 1:T){
          Lambda_t<-matrix(0,pall,R)
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
        D_z<-matrix(0,pall,pall)
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
        
        
        Phi<-matrix(0,pall,pall)
        for(t in 1:T){
          Lambda_t<-matrix(0,pall,R)
          for (i in 1:N) {
            Lambda_t<-Lambda_t+k(residual_tau[t,i]/hf)/hf*as.matrix(W[(i-1)*T+t,])%*%loading_tau[i,]
          }
          Lambda_t<-Lambda_t/N
          for (i in 1:N){
            Phi<-Phi+Lambda_t%*%solve(t(L_ols)%*%L_ols/N)%*%as.matrix(L_ols[i,])%*%X[(i-1)*T+t,]
          }
        }
        Phi<-Phi/(N*T)
        var_Mtheta<-matrix(0,T*N,pall)
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
        
      
        
        ###Computation of sum and variance (three terms)
        variance_tau<-matrix(0,pall,pall)
        variance_tau_robust<-matrix(0,pall,pall)
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
        
        
        t_size_c[,h_index]<-(beta_spjc[1:(pall)])/sqrt(diag(variance_tau)[1:(pall)])*sqrt(N*T)
        t_size_robust[,h_index]<-(beta_spjc[1:(pall)])/sqrt(diag(variance_tau_robust)[1:(pall)])*sqrt(N*T)
      }
      
      
      
      pvalue_c<-(1-pnorm(abs(t_size_c)))*2
      pvalue_robust<-(1-pnorm(abs(t_size_robust)))*2
      
      pvalue_c
      pvalue_robust

        b = matrix(0, nrow = pall-1, ncol = 4)
     

      
      b[, 1] = beta_spjc[2:length(beta_spjc)]
      b[, 2]  = sqrt(diag(variance_tau_robust)[2:(pall)])/sqrt(N*T)
      b[, 3] = t_size_robust[2:nrow(t_size_robust),1]
      b[, 4] =pvalue_robust[2:nrow(pvalue_robust),1]
      

      cnames = c( "Value", "Std. Error", "t value", "Pr(>|t|)")
      # if(Xid==1){
      #   rnames = c("intercept", "rho", "gdp_i", "gdp_j", "pop_i", "pop_j", "distance", "border", "language","europ", "ceep")
      # }
      if(Xid==1){
      if(dynamic=="T"){
          rnames = c( "rho", "gdp",  "pop","similarity", "distance", "border", "language","euro", "eec")
        } else{
          rnames = c( "gdp","pop","similarity", "distance", "border", "language","euro", "eec")
        }
      }

      # if(Xid==2){
      #   #rnames = c( "rho", "gdp",  "pop","similarity", "distance", "border", "language","Colony","euro", "eec")
      #   rnames = c( "dynamic", "gdp", "similarity", "distance", "border", "language","Colony","euro", "eec")
      # }
      # if(Xid==3|Xid==5){
      #   rnames = c( "rho", "gdp",  "pop","similarity", "distance", "border", "Colony","euro", "eec")
      # }
      # if(Xid==4){
      #   rnames = c( "rho", "gdp",  "pop","similarity", "rer","distance", "border", "language", "eec")
      # }
      # if(Xid==6){
      #   rnames = c( "rho", "gdp",  "pop","similarity", "distance", "border", "language","euro", "eec")
      # }
      # if(Xid==7){
      #   rnames = c( "rho", "gdp",  "pop","similarity", "distance", "border", "euro", "eec")
      # }
      # if(Xid==8){
      #   rnames = c( "rho", "gdp",  "pop","similarity", "distance", "border","language","colony", "euro", "eec")
      # }
      # 
      # if(Xid==9){
      #   rnames = c( "rho", "gdp",  "pop","similarity", "distance", "border", "colony","euro", "eec")
      # }
      # 
      # if(Xid==10){
      #   rnames = c( "rho", "gdp",  "pop","similarity","rer","rfl", "distance", "border", "language","euro", "eec")
      # }
      # 
      # if(Xid==11){
      #   rnames = c( "rho", "gdp",  "pop","similarity","rer","rfl", "distance", "border", "euro", "eec")
      # }
      # if(Xid==12){
      #   rnames = c( "rho", "gdp",  "pop","similarity","rfl", "distance", "border","language","euro", "eec")
      # }
      
      dimnames(b) = list(rnames, cnames)
      Netqrout<-list(b,TAU,R)
      
      Netallout[[q]] = Netqrout
    }
    
    # list.save(Netallout, 'Netallout.rds')
    filename<- paste(folder, '/Results/', 'Ydata_', YD, '_Xdata_', XD,"_FacN=",R, "_Tau=",taus[1],'_Correct=',correct,'_dynamic=',dynamic,'_(2)',sep ="")
    # 
    # saveRDS(Netallout, filename = sprintf("%s.rds", resultname))
    Netallout
    save(Netallout, file = paste(filename,".rds", sep=""))
    
    
    # Netallout<-readRDS("Netallout.rds")
    # str(Netallout)
    new<-list()
    for (li in 1:9){
      new[[li]]<-list("b"=Netallout[[li]][[1]],"tau"=Netallout[[li]][[2]])
    }
    # new
    # str(new)
    
    figname = paste(folder, '/Figure/', 'Ydata_', YD, '_Xdata_', XD,'_FacN=',R,'_Tau=',taus[1],'_Correct=',correct,'_dynamic=',dynamic,'_(2)','.pdf', sep = "")
    pdf(figname, paper="a4", height = 0, width = 0)
    Plotnetivrq(new)
    dev.off() 
  }
}



