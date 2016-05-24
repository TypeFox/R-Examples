##' Log-Likelihood for multiple ordinal right censored Multiple Ordinal Tobit (MOT) model.
##'
##' @title log-likelihood for mot model
##'
##' @param param parameter vector: (beta_0, beta_1, ... , beta_m, sigma).
##' @param xx design matrix of the model.
##' @param y observation vector.
##' @param tau threshold vector from tau_1 to tau_K.
##'
##' @return log-likelihood, vector with all observations.
##' @export 
##' @seealso \link[lmmot]{lmmot} 
##' @author Marvin Wright

motLogLik <- function(param,xx,y,tau) {
  x <- xx

  #sigma <- exp(param[length(param)])
	sigma <- param[length(param)]
  beta <- param[-length(param)]
	
	n <- length(y)
	K <- length(tau) 
	
	if (sigma < 0) {
		return(rep(-Inf, n))
	}
	
  yy <- t(x) %*% beta
  
  ll <- rep(NaN, n)
  
  # non censored data
  index <- y < tau[1]
  ll[index] <- dnorm((y[index] - yy[index])/sigma, log=TRUE) - log(sigma)
  
  # censored data, categories 1..K-1
  if (K > 1) {
    for (k in 1:(K-1)) {
      index <- (y >= tau[k] & y < tau[k+1])
      ll[index] <- log(pnorm((tau[k+1] - yy[index])/sigma) - 
        pnorm((tau[k] - yy[index])/sigma))
    }
  }
  
  # last category (K)
  index <- (y >= tau[K])
  ll[index] <- log(1 - pnorm((tau[K] - yy[index])/sigma))
  
  return(ll)
}