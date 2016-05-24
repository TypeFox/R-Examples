Max.Lam0 <- function(Y,J,B,in.lam0,alpha,thresh,maxiter,maxiter.B,
                     maxiter.T,nstep=10) {
  lam0 <- in.lam0
  est.B <- FLLat(Y,J,B,lam0*alpha,lam0*(1-alpha),thresh,maxiter,maxiter.B,
                 maxiter.T)$Beta
  check <- sum((scale(est.B,scale=F))^2)
  min.max <- c(0,Inf)
  for (i in 1:nstep) {
    if (check<=10^(-10)) {
      min.max[2] <- lam0
      lam0 <- mean(min.max)
      est.B <- FLLat(Y,J,B,lam0*alpha,lam0*(1-alpha),thresh,maxiter,
                     maxiter.B,maxiter.T)$Beta
      check <- sum((scale(est.B,scale=F))^2)
    } else {
      min.max[1] <- lam0
      if (min.max[2]==Inf) {
        lam0 <- 2*lam0
      } else {
        lam0 <- mean(min.max)
      }
      est.B <- FLLat(Y,J,B,lam0*alpha,lam0*(1-alpha),thresh,maxiter,
                     maxiter.B,maxiter.T)$Beta
      check <- sum((scale(est.B,scale=F))^2)
    }
  }
  if (check<=10^(-10)) {
    return(lam0)
  } else {
    return(min.max[2])
  }
}

