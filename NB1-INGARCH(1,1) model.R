rm(list = ls())
gc()#
#Log-likelihood of Zhu's NB-INGARCH(1,1) model:
llzhunbingarch11 <- function(par,n=3,data){
  #par is vector (beta0,beta1,alpha1,n,m1)
  beta0 <- par[1]
  beta1 <- par[2]
  alpha1 <- par[3]
  cmean <- mean(data) #Initialization
  T <- length(data)
  value <- 0 #conditional likelihood
  for(t in c(2:T)) {
    cmean <- beta0+alpha1*data[t-1]+beta1*cmean
    p <- 1/(1+cmean/n)
    value <- value-log(dnbinom(data[t],n,p))
  }
  value
}
ui = matrix(0, nrow = 4, ncol = 3)
ui[1, 1] = 1
ui[2, 2] = 1
ui[3, 3] = 1
ui[4, 2] = -1
ui[4, 3] = -1
ci = rep(0, 4)
ci[4] = -1
data <- read.table(file="C:\\Users\\Administrator\\Desktop\\Downloads.txt",
header=F, encoding="UTF-8")
Tlen <- length(data$V1)
mean(data$V1)
var(data)
rho1 <- acf(data, plot=FALSE)[[1]][2]
rho1 
rho2 <- acf(data, plot=FALSE)[[1]][3]
rho2 
barX <- mean(data$V1)
delta1mm <- rho2/rho1
beta0mm <- barX*(1-delta1mm)
alpha1mm <- (1-delta1mm^2 - sqrt(1-delta1mm^2) * sqrt(1-delta1mm^2+4*rho1*(delta1mm-rho1)) )/(-2*(delta1mm-rho1))
beta1mm <- delta1mm-alpha1mm
theta_init=c(beta0mm,alpha1mm,beta1mm)
estmlnbz <- constrOptim(theta=theta_init,f=llzhunbingarch11,data=data$V1,
ui=ui,ci=ci,outer.iterations = 100,outer.eps = 1e-05,method = "Nelder-Mead") 
beta0mlnbz <- estmlnbz$par[[1]]
alpha1mlnbz <- estmlnbz$par[[2]]
beta1mlnbz <- estmlnbz$par[[3]]
ofiestnbz <- estmlnbz$hessian 
neglmaxnbz <- estmlnbz$value
#Estimates:
c(beta0mlnbz,alpha1mlnbz,beta1mlnbz)
#AIC and BIC:
AICnbz <- 2*neglmaxnbz+2*4
BICnbz <- 2*neglmaxnbz+log(Tlen)*4
c(neglmaxnbz, AICnbz, BICnbz)
 
