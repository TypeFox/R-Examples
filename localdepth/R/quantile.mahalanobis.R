#############################################################
#
#	quantile.mahalanobis function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: July, 19, 2011
#	Version: 0.3-1
#
#	Copyright (C) 2011 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

quantile.mahalanobis <- function(x, probs, nsamp='all', all=FALSE, covariance=NULL, ...) {
  x <- data.matrix(x)
  if (ncol(x) < 2) stop('At least two columns must be supplied')
  nx <- nrow(x)
  if (is.null(covariance))
    covariance <- cov(x)
  if (is.numeric(nsamp) && nsamp <= 0) stop("the argument 'nsamp' must be positive")
  if (is.numeric(nsamp)) {
    index <- sample(1:nx, size=nsamp, replace =TRUE)
    nx <- nsamp
    x <- x[index,]
  }  
  if (abs(det(covariance)/2) < .Machine$double.eps) stop('The covariance matrix seems to be singular')
  Sinv <- solve(covariance)
  mah <- rep(0, (nx*(nx-1))/2) 
  k <- 0
  for(i in 1:(nx-1)) {
    for(j in (i+1):nx) {
      k <- k+1
      temp <- x[i,]-x[j,]
      mah[k] <- sqrt(temp%*%Sinv%*%temp)
    }
  }
  res <- quantile(mah, probs, ...)
  if (all)
    res <- list(quantile=res, stats=mah, call=match.call())
  class(res) <- 'quantile.localdepth'
  return(res)
}
