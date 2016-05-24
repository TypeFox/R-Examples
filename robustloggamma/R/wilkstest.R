#############################################################
#	All functions in this file are copyrighted
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: July, 10, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

loggammarob.wilks <- wWilksTest <- function(x, thetainit=NULL, method="L-BFGS-B", lower=c(-Inf, 0.0001), upper=c(Inf, Inf), ...) {
  y <- x$data
  if (is.null(x$weights)) 
    weights <- rep(1, length(y))
  else
    weights <- x$weights
  
  minusloglikH0 <- function(theta, x, weights) {
    if (is.null(dim(theta)))
      res <- drop(dloggamma(x=x, mu=theta[1], sigma=theta[2], lambda=theta[2], log=TRUE)%*%weights)
    else {
      res <- rep(NA, NROW(theta))    
      for (i in 1:NROW(theta)) {
        res[i] <- drop(dloggamma(x=x, mu=theta[i,1], sigma=theta[i,2], lambda=theta[i,2], log=TRUE)%*%weights)
      }
    }
    res <- -res
    return(res)
  }
  if (is.null(thetainit))
    thetainit <- c(x$mu, x$sigma)

  theta0 <- optim(par=thetainit, fn=minusloglikH0, x=y, weights=weights, method=method, lower=lower, upper=upper, ...)$par
  
  statistic <- drop(2*(dloggamma(x=y, mu=x$mu, sigma=x$sigma, lambda=x$lambda, log=TRUE) - dloggamma(x=y, mu=theta0[1], sigma=theta0[2], lambda=theta0[2], log=TRUE))%*%weights)
  result <- list()
  names(statistic) <- "wilks"
  result$statistic <- statistic
  parameter <- 1
  names(parameter) <- "df"
  result$parameter <- parameter
  result$p.value <- pchisq(statistic, df=1, lower.tail=FALSE)
  result$conf.int <- NULL
  names(theta0) <- c("mu", "sigma and lambda")
  result$estimate <- theta0
  result$null.value <- "true lambda"
  names(result$null.value) <- "sigma"
  result$alternative <- "two.sided"
  result$method <- "weighted Wilks test"
  result$data.name <- paste(deparse(substitute(x$data)))
  result$mu0 <- theta0[1]
  result$sigma0 <- theta0[2]
  result$lambda0 <- theta0[2]
  class(result) <- "htest"
  return(result)
}
