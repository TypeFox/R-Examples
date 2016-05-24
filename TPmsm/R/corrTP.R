corrTP <- function(dist, corr, dist.par) {
	if ( missing(dist) ) stop("Argument 'dist' is missing, with no default")
	if ( missing(corr) ) stop("Argument 'corr' is missing, with no default")
	if ( missing(dist.par) ) stop("Argument 'dist.par' is missing, with no default")
	if ( !( dist %in% c("weibull", "exponential") ) ) stop("Argument 'dist' must be one of 'weibull' or 'exponential'")
	if (dist == "exponential") {
		if (corr < -1 | corr > 1) stop("Argument 'corr' with dist='exponential' must be greater or equal to -1 and lower or equal to 1")
		if (length(dist.par) != 2) stop("Argument 'dist.par' with 'dist=exponential' must be a vector with lenght 2")
		if (dist.par[1] <= 0 | dist.par[2] <= 0) stop("Argument 'dist.par' must be greater than 0")
		corrxy <- (1/4)*corr
		return(corrxy)
	}
	if (dist == "weibull") {
		if (corr <= 0 | corr > 1) stop("Argument 'corr' with 'dist=weibull' must be greater than 0 and lower or equal to 1")
		if (length(dist.par) != 4) stop("Argument 'dist.par' with 'dist=weibull' must be a vector with lenght 4")
		if (dist.par[1] <= 0 | dist.par[2] <= 0 | dist.par[3] <= 0 | dist.par[4] <= 0) stop("Argument 'dist.par' must be greater than 0")
		a <- gamma( (corr/dist.par[1])+1 )
		b <- gamma( (corr/dist.par[3])+1 )
		c <- gamma( (1/dist.par[1])+(1/dist.par[3])+1 )
		d <- gamma( (1/dist.par[1])+1 )
		e <- gamma( (1/dist.par[3])+1 )
		f <- gamma( (corr/dist.par[1])+(corr/dist.par[3])+1 )
		cov1 <- dist.par[2]*dist.par[4]*(a*b*c-d*e*f)/f
		varx <- dist.par[2]^2*(gamma( (2/dist.par[1])+1 )-d^2)
		vary <- dist.par[4]^2*(gamma( (2/dist.par[3])+1 )-e^2)
		corrxy <- cov1/sqrt(varx*vary)
		return(corrxy)
	}
}
