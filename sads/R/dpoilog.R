dpoilog <- function(x, mu, sig, log = FALSE){
	#### FIX: poilog::dpoilog throws an error if an invalid parameter is entered
	#### so we have to circumvent the error here:
	if (length(mu) > 1 | length(sig) > 1) stop ("Vectorization of parameters not implemented")
	if (any(!is.wholenumber(x))) warning("non integer values in x")
	to.zero <- ! is.wholenumber(x) | x < 0
	to.NaN <- NULL
	if (!is.finite(mu)) to.NaN <- 1:length(x)
	if (!is.finite(sig) | sig <= 0) to.NaN <- 1:length(x)
	x[! is.wholenumber(x) | x < 0] <- 0
	mu[ !is.finite(mu) ] <- 1
	sig[ !is.finite(sig) | sig <= 0] <- 1
	y <- poilog::dpoilog(x, mu, sig)
	y[to.NaN] <- NaN
	y[to.zero] <- 0
	if (any(is.nan(y))) warning ("NaNs produced")
	if (log) return(log(y))
	else return(y)
}
