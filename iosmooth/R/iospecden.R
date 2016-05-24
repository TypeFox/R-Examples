iospecden <-
function(x, l, kernel = c("Trap", "Rect", "SupSm"), x.points = seq(-pi, pi, len=200)) {
	# x is the time series
	# l is the bandwidth
	kernel <- match.arg(kernel, c("Trap", "Rect", "SupSm"))

	kappa <- switch(kernel, Trap=kappaTrap, Rect=kappaRect, SupSm=kappaInDf)
	
	if(missing(l)) l <- bwadap.ts(x)
	xacf <- as.vector(acf(x, lag.max=max((2*l-1), 0), type="cov", plot=F)$acf)
		
	Fones <- function(s) {
		# spectral density at a single s
		# Formula below takes real part at lags 1 to 2l - 1, doubles it
		# then adds a single autocovariance at lag 0  
		(2*sum(kappa(1:(2*l-1), l) * xacf[-1] * cos(s * 1:(2*l-1))) + xacf[1])/(2*pi)   
	}

	f <- Vectorize(Fones, "s")
	if( is.null(x.points) ) {
		return(f) #Returns the vectorized spectral density function
	}
	return(list(x=x.points, y=f(x.points)))
}
