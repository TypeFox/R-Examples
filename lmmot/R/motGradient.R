##' Gradient of log-Likelihood for right censored Multiple Ordinal Tobit (MOT) model.
##'
##' @title Gradient of log-Likelihood for mot model
##'
##' @param param parameter vector: (beta_0, beta_1, ... , beta_m, sigma).
##' @param xx design matrix of the model.
##' @param y observation vector.
##' @param tau threshold vector from tau_1 to tau_K.
##'
##' @return gradient of log-likelihood, vector with all observations.
##' @seealso \link[lmmot]{lmmot} 
##' @author Marvin Wright

motGradient <- function(param,xx,y,tau) {
  x <- xx

  #sigma <- exp(param[length(param)])
	sigma <- param[length(param)]
	beta <- param[-length(param)]
  
  n <- length(y)
  m <- length(param)
  K <- length(tau)
  
  yy <- t(x) %*% beta
  
  grad <- matrix(rep(NaN, n*m), n, m)

  # non censored data
  index <- y < tau[1]
	if (sum(index) > 0) {
    dbeta <- 1/sigma^2 * t(x[,index]) * (y[index] - yy[index]) 
    dsigma <- 1/sigma^3 * (y[index] - yy[index])^2 - 1/sigma
    grad[index,] <- cbind(dbeta, dsigma)
	}
 
  # censored data, categories 1..K-1
  if (K > 1) {
    for (k in 1:(K-1)) {
      index <- (y >= tau[k] & y < tau[k+1])
      if (sum(index) > 0) {
        dbeta <- t(x[,index])/sigma * 
          ( dnorm((tau[k]-yy[index])/sigma) -
            dnorm((tau[k+1]-yy[index])/sigma) ) /
          ( pnorm((tau[k+1]-yy[index])/sigma) -
            pnorm((tau[k]-yy[index])/sigma) )
        dsigma <- 1/sigma^2 * 
          ( (tau[k]-yy[index]) * dnorm((tau[k]-yy[index])/sigma) -
            (tau[k+1]-yy[index]) * dnorm((tau[k+1]-yy[index])/sigma) ) / 
          ( pnorm((tau[k+1]-yy[index])/sigma) - 
            pnorm((tau[k]-yy[index])/sigma) )
        grad[index,] <- cbind(dbeta, dsigma)
      }
    }
  }
  
  # last category (K)
  index <- (y >= tau[K])
	if (sum(index) > 0) {
    dbeta <- t(x[,index])/sigma * 
      dnorm((tau[K]-yy[index])/sigma) /
      (1 - pnorm((tau[K]-yy[index])/sigma))
    dsigma <- 1/sigma^2 * 
      (tau[K]-yy[index]) * 
      dnorm((tau[K]-yy[index])/sigma)  / 
      (1 - pnorm((tau[K]-yy[index])/sigma))
    grad[index,] <- cbind(dbeta, dsigma)
	}

  return(grad)
  
}