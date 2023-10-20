m(list = ls())
gc()
library("MASS")
#infromation matri bases on the poisson distribution
information.poisson <- function(theta,data){
  n <- length(data)
  lambda=rep(NA,n)
  lambda[1]=mean(data)
  first=rep(NA,n)
  first[1]=1
  second=rep(NA,n)
  second[1]=1
  third=rep(NA,n)
  third[1]=1
  Information=matrix(9,nrow = 3,ncol = 3)
  s1=rep(NA,n)
  s2=rep(NA,n)
  s3=rep(NA,n)
  for (t in c(2:n)){
    lambda[t]=theta[1]+theta[2]*lambda[t-1]+theta[3]*data[t-1]
    first[t]=(1+theta[2]*first[t-1])
    second[t]=(lambda[t-1]+theta[2]*second[t-1])
    third[t]=(data[t-1]+theta[2]*third[t-1])
    s1[t]=first[t]
    s2[t]=second[t]
    s3[t]=third[t]
    var.comp=(1/sqrt(lambda[t]))*c(s1[t],s2[t],s3[t])
    Information=Information+var.comp%*%t(var.comp)
  }
  return(Information)
}
#information based on negative binmial distribution
information.negbin <- function(theta,data,nu){
  n <- length(data)
  lambda=rep(NA,n)
  lambda[1]=mean(data)
  first=rep(NA,n)
  first[1]=1
  second=rep(NA,n)
  second[1]=1
  third=rep(NA,n)
  third[1]=1
  Information=matrix(9,nrow = 3,ncol = 3)
  s1=rep(NA,n)
  s2=rep(NA,n)
  s3=rep(NA,n)
 for (t in c(2:n)){
    lambda[t]=theta[1]+theta[2]*lambda[t-1]+theta[3]*data[t-1]
    first[t]=(1+theta[2]*first[t-1])
    second[t]=(lambda[t-1]+theta[2]*second[t-1])
    third[t]=(data[t-1]+theta[2]*third[t-1])
    s1[t]=first[t]
    s2[t]=second[t]
    s3[t]=third[t]
    var.comp=(sqrt(1/lambda[t]+1/nu))*c(s1[t],s2[t],s3[t])
    Information=Information+var.comp%*%t(var.comp)
  }
  return(Information)
}
#Log-likelihood of Poisson INGARCH(1,1) model:
llpingarch11 <- function(par,data){
  #par is vector (beta0,beta1,alpha1,m1)
  beta0 <- par[1]
  beta1 <- par[2]
  alpha1 <- par[3]
  cmean <- mean(data) #Initialization of conditional mean
  T <- length(data)
  value <- 0 #conditional likelihood
  for(t in c(2:T)) {
    cmean <- beta0+alpha1*data[t-1]+beta1*cmean
    value <- value-log(dpois(data[t],cmean))
  }
  value
}
data <- read.table(file="C:\\Users\\Administrator\\Desktop\\Downloads.txt",
header=F, encoding="UTF-8")
Tlen <- length(data$V1)
rho1 <- acf(data, plot=FALSE)[[1]][2]
rho1 
rho2 <- acf(data, plot=FALSE)[[1]][3]
rho2 
barX <- mean(data$V1)
ui = matrix(0, nrow = 4, ncol = 3)
ui[1, 1] = 1
ui[2, 2] = 1
ui[3, 3] = 1
ui[4, 2] = -1
ui[4, 3] = -1
ci = rep(0, 4)
ci[4] = -1
barX <- mean(data$V1)
delta1mm <- rho2/rho1
beta0mm <- barX*(1-delta1mm)
alpha1mm <-(1-delta1mm^2 - sqrt(1-delta1mm^2) * sqrt(1-delta1mm^2+4*rho1*(delta1mm-rho1)) )/(-2*(delta1mm-rho1))
beta1mm <- delta1mm-alpha1mm
theta_init=c(beta0mm,alpha1mm,beta1mm)
estmlp <- constrOptim(theta=theta_init,f=llpingarch11,data=data$V1,
ui=ui,ci=ci,outer.iterations = 100,outer.eps = 1e-05,method = "Nelder-Mead") 
beta0mlp <- estmlp$par[1]
alpha1mlp <- estmlp$par[2]
beta1mlp <- estmlp$par[3]
neglmaxp <- estmlp$value
#Estimates:
round(c(beta0mlp,alpha1mlp,beta1mlp),5)
#AIC and BIC:
AICp <- 2*neglmaxp+2*5
BICp <- 2*neglmaxp+log(Tlen)*5
c(neglmaxp, AICp, BICp)
hat.d <- beta0mlp
hat.a <- alpha1mlp
hat.b <- beta1mlp
lambda <- rep(NA,Tlen)
lambda[1] <- mean(data$V1)
for (t in 2:267){
  lambda[t] <- hat.d+hat.b*lambda[t-1]+hat.a*data$V1[t-1]
}
#Estimation \nu
est.nu1 <- (mean(((data$V1-lambda)^2-lambda)/(lambda^2)))^(-1)
est.nu1 
est.nu2 <- theta.mm(data$V1,lambda,Tlen-3)
est.nu2 
