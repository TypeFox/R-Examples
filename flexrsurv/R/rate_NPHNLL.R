rate_NPHNLL <- function(gamma, beta0, beta, alpha, T0, T, X, Z, ntd, nnll, nsbtd, nsbnll, time.spline, z.spline ){
  # compute the rate of mahaboubi model
  return(exp(T0 %*% gamma) + X %*% beta0 + T %*% beta %*% t(Z %*% alpha)) 
}

