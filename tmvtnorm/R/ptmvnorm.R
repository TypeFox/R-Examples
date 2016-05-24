
# Verteilungsfunktion der truncated multivariate normal distribution
#
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
ptmvnorm <- function(lowerx, upperx, mean=rep(0, length(lowerx)), sigma, lower = rep(-Inf, length = length(mean)), upper = rep( Inf, length = length(mean)), maxpts = 25000, abseps = 0.001, releps = 0)
{
  # check of standard tmvtnorm arguments
  cargs <- checkTmvArgs(mean, sigma, lower, upper)
  mean  <- cargs$mean
  sigma <- cargs$sigma
  lower <- cargs$lower
  upper <- cargs$upper
  
  # check of additional arguments lowerx and upperx
  if (is.null(lowerx) || any(is.na(lowerx))) 
    stop(sQuote("lowerx"), " not specified or contains NA")
  if (is.null(upperx) || any(is.na(upperx))) 
    stop(sQuote("upperx"), " not specified or contains NA")  
  if (!is.numeric(lowerx) || !is.vector(lowerx)) 
    stop(sQuote("lowerx"), " is not a numeric vector")
  if (!is.numeric(upperx) || !is.vector(upperx)) 
    stop(sQuote("upperx"), " is not a numeric vector")  
  if (length(lowerx) != length(lower) || length(lower) != length(upperx))
	stop("lowerx an upperx must have the same length as lower and upper!")  
  if (any(lowerx>=upperx))
	stop("lowerx must be smaller than or equal to upperx (lowerx<=upperx)")
    
  # Aufpassen: 
  # Wir müssen garantieren, dass nur innerhalb des Support-Bereichs lower <= x <= upper integriert wird. Sonst kann Ergebnis >= 1 rauskommen.
  # Wenn einzelne Komponenten von lowerx <= lower sind, dann von der Untergrenze lower integrieren. Analog für upperx >= upper
  f <- pmvnorm(lower=pmax(lowerx, lower), upper=pmin(upperx, upper), mean=mean, sigma=sigma, maxpts = maxpts, abseps = abseps, releps = releps) / 
	   pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma, maxpts = maxpts, abseps = abseps, releps = releps)
  return(f)
}
