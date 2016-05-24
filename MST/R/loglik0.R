loglik0 <-
function(para, dat){ 
  beta0 <- para[1]; v <- para[2]; Delta.i <- dat$Delta.i; X.i <- dat$X.i
  n <- length(X.i)
  A.i <- 1/v + Delta.i
  B.i <- 1/v + exp(beta0)*X.i
  options(warn=-1)
  loglik0 <- beta0*sum(Delta.i) + sum(log(gamma(A.i))) - n*log(gamma(1/v)) - n/v*log(v) - sum(A.i*log(B.i))
  options(warn=0)
  return(-loglik0)
}
