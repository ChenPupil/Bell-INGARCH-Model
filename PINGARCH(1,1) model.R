rm(list = ls())
gc()
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
delta1mm <- rho2/rho1
beta0mm <- barX*(1-delta1mm)
alpha1mm <- (1-delta1mm^2 - sqrt(1-delta1mm^2) * sqrt(1-delta1mm^2+4*rho1*(delta1mm-rho1)) )/(-2*(delta1mm-rho1))
beta1mm <- delta1mm-alpha1mm
ui = matrix(0, nrow = 4, ncol = 3)
ui[1, 1] = 1
ui[2, 2] = 1
ui[3, 3] = 1
ui[4, 2] = -1
ui[4, 3] = -1
ci = rep(0, 4)
ci[4] = -1
estmlp <- constrOptim(theta=c(beta0mm,alpha1mm,beta1mm),f=llpingarch11,
data=data$V1,ui=ui,ci=ci,outer.iterations = 100,outer.eps = 1e-05,
method = "Nelder-Mead") 
beta0mlp <- estmlp$par[[1]]
alpha1mlp <- estmlp$par[[2]]
beta1mlp <- estmlp$par[[3]]
ofiestp <- estmlp$hessian 
neglmaxp <- estmlp$value
#Estimates:
round(c(beta0mlp,alpha1mlp,beta1mlp),3)
#AIC and BIC:
AICp <- 2*neglmaxp+2*3
BICp <- 2*neglmaxp+log(Tlen)*3
c(neglmaxp, AICp, BICp)

