dpoig <- function(x, frac, rate, shape, log=FALSE) {
	frac[ !is.finite(frac) | frac <= 0 ] <- NaN
	rate[ !is.finite(rate) | rate <= 0 ] <- NaN
	shape[ !is.finite(shape) | shape <= 0 ] <- NaN
	b <- x*log(frac)+shape*log(rate)+lgamma(x+shape)
	c <- lfactorial(x)+lgamma(shape)+(x+shape)*log(frac+rate)
	vals<-b-c
	if (any(is.nan(vals))) warning ("NaNs produced")
	if (any(!is.wholenumber(x))) warning("non integer values in x")
	vals[ ! is.wholenumber(x) | x < 0] <- -Inf
	if(log)vals else exp(vals)
}
