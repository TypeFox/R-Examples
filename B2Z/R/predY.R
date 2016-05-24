predY <-
function(parms,n)
  {
  sigma <- matrix(c(parms[1],parms[2],parms[2],parms[3]), 2, 2)
  m <- matrix(parms[-c(1,2,3)],n,2)
  ev <- eigen(sigma, symmetric = TRUE)
  retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
  retval <- m + matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
  return(as.numeric(retval))
  }

