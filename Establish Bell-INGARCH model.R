#################################
ts <- function(n,d,a,b) {
  x = rep(0,n)
  mu = rep(0,n)
  mu[1] = 1 # intial value
  x[1] = rbell(1,mu[1])
    for (t in c(2:n)) {
      mu[t] = d+a*x[t-1]+b*mu[t-1]
      x[t] = rbell(1,mu[t])
    }
    return(x)
  }
