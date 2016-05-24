#############################################################
#	All functions in this file are copyrighted
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

mle.gamma <- function(x, tol=1e-8) {
  n <- length(x)
  shape1 <- mean(x)^2/(mean(x^2)-mean(x)^2)
  sumlogx <- sum(log(x))
  logmeanx <- log(mean(x))

  g <- function(x) sumlogx + n*(log(x) - logmeanx) - n*digamma(x)
  gprime <- function(x) n/x - n*trigamma(x) 
  nr <- function(x) x - g(x)/gprime(x)
  diffshape <-  tol + 1
  while (abs(diffshape) > tol) {
    shape2 <- nr(shape1)
    diffshape <- shape2 - shape1
    shape1 <- shape2
  }
  rate <- shape1/mean(x)
  ans <- list(shape=shape2, rate=rate, scale=1/rate)
  return(ans)
}
