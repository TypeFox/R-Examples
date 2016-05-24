derivs <- function(tnow, y, bwtlist) {
	#  Sets up a linear differential equation of order m
	#  as a system of m first order equations
	#  Arguments:
	#  TNOW    ... A vector of values of the independent variable t
	#  Y       ... A matrix of m values of Y corresponding to TNOW
	#  BWTLIST ... A functional data object containing coefficient functions
	#  Returns:
	#  DY      ... A matrix of derivative values corresponding to Y
	
	#  Last modified:  26 October 2005
	
  m  <- length(bwtlist);
  wmat <- matrix(0, m, m)
  wmat[1:(m-1),2:m] <- diag(rep(1,m-1))
  for (j in 1:m) {
	   bfdParj <- bwtlist[[j]]
	   wj      <- eval.fd(tnow, bfdParj$fd)
	   wmat[m,j] <- -wj
  }
  dy <- wmat %*% y
  return(dy)
}
