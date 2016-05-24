dpoix <- function(x, frac, rate, log=FALSE) {
	frac[ !is.finite(frac) | frac <= 0 ] <- NaN
	rate[ !is.finite(rate) | rate <= 0 ] <- NaN
	b <- x*log(frac)
	m <- log(rate)
	n <- (x+1)*log(rate+frac)
	vals <- b+m-n
	if (any(is.nan(vals))) warning ("NaNs produced")
	if (any(!is.wholenumber(x))) warning("non integer values in x")
	vals[ ! is.wholenumber(x) | x < 0] <- -Inf
	if(log) vals else exp(vals)
}


