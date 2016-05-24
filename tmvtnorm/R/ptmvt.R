
# Verteilungsfunktion der truncated multivariate t distribution
#
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
ptmvt <- function(
    lowerx, 
    upperx, 
    mean=rep(0, length(lowerx)),
    sigma, df = 1, 
		lower = rep(-Inf, length = length(mean)), 
    upper = rep( Inf, length = length(mean)), 
		maxpts = 25000, abseps = 0.001, releps = 0)
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
	f <- pmvt(lower=pmax(lowerx, lower), upper=pmin(upperx, upper), delta=mean, sigma=sigma, df=df, maxpts = maxpts, abseps = abseps, releps = releps, type="shifted") / 
			pmvt(lower=lower, upper=upper, delta=mean, sigma=sigma, df=df, maxpts = maxpts, abseps = abseps, releps = releps, type="shifted")
	return(f)
}