###### d,a,b>0,a+b<1 #####
ui = matrix(0, nrow = 4, ncol = 3)
ui[1, 1] = 1
ui[2, 2] = 1
ui[3, 3] = 1
ui[4, 2] = -1
ui[4, 3] = -1
ci = rep(0, 4)
ci[4] = -1

#####  CMLE #####
cmle <- function(d,a,b,n,reps){
est1 <- matrix(NA,nrow=reps,ncol=3)
MADE <- matrix(NA,nrow=reps,ncol=3)
  for (i in c(1:reps)){
    #Estimation
    daten <- ts(n, d, a, b)
    d1 <- runif(1,d-0.001,d+0.001)
    a1 <- runif(1,a-0.001,a+0.001)
    b1 <- runif(1,b-0.001,b+0.001)
    theta_init=c(d1,a1,b1)
    vp <- constrOptim(theta=theta_init,f=GBELLloglike,data=daten,ui=ui,ci=ci,
    outer.iterations = 100,outer.eps = 1e-05,method = "Nelder-Mead")$par                       
    dmad <- abs(vp[1]-d);amad <- abs(vp[2]-a);bmad <- abs(vp[3]-b);
    MADE[i,] <- c(dmad,amad,bmad)
    est1[i,] <- vp
  }
  estCMLE <- apply(est1,2,mean)   #
  biasCMLE <- estCMLE-c(d,a,b)    #
  madeCMLE <- apply(MADE,2,mean)  #
  varCMLE <- apply(est1,2,var)    #
  mseCMLE <- varCMLE+biasCMLE^2   #
  cat("N=", n,"\n")
  cat("d=", d,"\n")
  cat("a=", a,"\n")
  cat("b=", b,"\n")
  cat("#---------------------------MEAN------------------------#","\n")
  print.est <- rbind(estCMLE)
  rownames(print.est) <- c("CMLE")
  colnames(print.est) <- c("d","a","b")
  print(round(print.est,4))
  cat("#---------------------------BIAS------------------------#","\n")
  print.est <- rbind(biasCMLE)
  rownames(print.est) <- c("CMLE")
  colnames(print.est) <- c("d","a","b")
  print(round(print.est,4))
  cat("#---------------------------MADE------------------------#","\n")
  print.est <- rbind(madeCMLE)
  rownames(print.est) <- c("CMLE")
  colnames(print.est) <- c("d","a","b")
  print(round(print.est,4))
  cat("#---------------------------MSE------------------------#","\n")
  print.est <- rbind(mseCMLE)
  rownames(print.est) <- c("CMLE")
  colnames(print.est) <- c("d","a","b")
  print(round(print.est,4))
}
