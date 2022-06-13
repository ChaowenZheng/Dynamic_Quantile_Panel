
####Code for dynamic quantile panel regression empircial analysis


rm(list=ls())


library(quantreg,lib.loc = "H:/rlibs/4.0.2")
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

#setwd("//userfs/cz1113/w2k/Desktop/Chapter2")
#setwd("\\Ecopc1055\c$\Users/cz1113/Google Drive/PhD/Potential Topic/Quantile/Dynamic Quantile Panel Data with Interactive Effects")

#===================================#
##          Define Functions       ##
#===================================#

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
rq.fit.PanelIE <- function(X_final,AY,tau=TAU,R,iter=200,cri_V=0.00001,c1,power1,smooth){
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
  
  ###Second, substitute the factors and run quantile regression
  
  IN<-diag(1,n,n)
  IT<-diag(1,t,t)
  IF<-kronecker(IN, F)
  
  #Alternative way (This seems to be the best way)
  #Instead of defactor first, we directly write the model using dummy variable
  #We therefore could obtain quantile dependent coefficient for beta
  #quantile and individual specific coefficient for factor loadings
  y_final<- vec(AY[2:nrow(AY),])
  X_final<-cbind(vec(AY[1:(nrow(AY)-1),]), X_final,IF) 
  
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

#===================================#
##               Design            ##
#===================================#
#setwd("//Ecopc1055/c$/Users/cz1113/Google Drive/PhD/chapter2_empirical")
setwd("C:/Users/86188/Google ÔÆ¶ËÓ²ÅÌ/PhD/chapter2_empirical")
#data<-read.csv("eu2020.csv",header = T)
data<-read.csv("oecd2020.csv",header = T)
country<-c(unique(data$codeh),"USA")
c_ind<-0
res_all<-matrix(0,20,6)
for (coun in (country)){
  c_ind<-c_ind+1
data<-read.csv("oecd2020.csv",header = T)
data<-data[which(data$codeh==coun|data$codef==coun), ]
# Number of individuals
N<-length(unique(data$pair_n))
IN<-diag(1,N,N)
# data$gdppc<-log(data$gdppcf+data$gdppch)
# data$pop<-log(data$popf+data$poph)
# write.csv(data,"eu2020.csv")
#Time periods
T<-max(data$year)-min(data$year)
#IT<-diag(1,T,T)
#Trend
trend<-rep(1:(T+1),each=N)
#Specifying y
#Y<-data$trade
Y1<-log(data$rexport)#+log(data$rimport)
Y2<-log(data$rimport)
#Y<-log(Y)
#Y<-data$rimport
#Specifying x
X<-cbind(log(data$gdpf/100)+log(data$gdph/100),data$sim,data$rer,data$rfl,data$ceep)
X<-cbind(data$gdp,data$sim,data$rer,data$rfl,data$ceep)
X<-cbind(log(data$gdpf),log(data$gdph),log(data$popf),log(data$poph),log(data$distance),data$border,data$ceep)
X<-cbind(log(data$gdpf),log(data$gdph),log(data$popf),log(data$poph),log(data$distance))
#X<-cbind(log(data$gdpf),log(data$gdph),log(data$popf),log(data$poph),data$border,data$language,log(data$distcap),data$sim,data$rer,data$rfl,data$europ,data$ceep)
#X<-cbind(data$gdp,data$sim,data$rer,data$rfl,data$ceep,trend)
#X<-cbind(data$gdp,data$sim,data$rer,data$rfl,data$ceep,data$pop)
#X<-cbind(data$gdp,data$rer,data$rfl,data$ceep)
dim(X)
p<-ncol(X)
matplot(X, type = "l")


Y1_wide<-t(matrix(Y1,N,T+1,))
Y1_wide<-Y1_wide[1:(T+1),]

Y2_wide<-t(matrix(Y2,N,T+1,))
Y2_wide<-Y2_wide[1:(T+1),]

#Y_wide<-t(matrix(fit$residuals,N,T+1,))
#FOR GROWTH
#Y_wide<-diff(log(Y_wide))
# matplot(Y_wide[,1:5], type = "l")
#for level data
X_wide<-array(0,dim=c(T+1,N,p))
for (i in 1:p){
  x_tem<-t(matrix(X[,i],N,T+1))
  X_wide[,,i]<-x_tem
}
dim(X_wide)
#dim(Y_wide)
# #for growth data
# X_wide<-array(0,dim=c(T,N,p))
# for (i in 1:p){
#   if(i ==1){
#     x_tem<-diff(t(matrix(X[,i],N,T+1)))
#   } else{
#     x_tem<-t(matrix(X[,i],N,T+1))
#     x_tem<-x_tem[2:(T+1),]
#   }
#   X_wide[,,i]<-x_tem
# }
X_wide<-X_wide[2:(T+1),,]
# # for diff
# X_wide<-X_wide[2:T,,]
# dim(X_wide)
# #minus lag or diff
# T<-T-1


X_final<-matrix(0,N*(T),p)
for (i in 1:p){
  x_tem<-matrix(X_wide[,,i],N*(T),1)
  X_final[,i]<-x_tem
}
matplot(X_final, type = "l")
#X_final

iter<-200 # maximum iteration for IPC
cri_V<-10^-7 # Convergence criterion
TAU<-0.5
R<-2

Y_wide<-cbind(Y1_wide,Y2_wide)
X_final<-rbind(X_final,X_final)
# dim(Y_wide)
# dim(X_final)
#===================================#
#Estimation 
#===================================#

########   Our Proposed estimator

###  First, run IPC to obtain estimator for factors
#In Galvao(2016), the bandwidth is choosen based on the variance of residual of non-smoothed quantile regression
##For bandwidth in the smoothed regression
c1<-2
power1<--1/10
fit_all<-rq.fit.PanelIE(X_final,Y_wide,TAU,R,iter,cri_V,c1,power1, smooth="F")
fit<-fit_all$fit
#b_smooth<-t(as.matrix(fit_all$sfit[2:4]))
b <- fit$coefficients[1:(p+2)]
b
res_all[c_ind,]<-b[2:7]
print(c_ind)
print(dim(Y_wide))
}
rownames(res_all)<-country
res_all