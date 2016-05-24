qpoig <- function(p, frac, rate, shape, S=30, lower.tail=TRUE, log.p=FALSE) {
	if (log.p) p <- exp(p)
	if(!lower.tail) p <- 1 - p
	q <- function(p) suppressWarnings(qfinder("poig", p, list(frac=frac, rate=rate, shape=shape)))
	y <- sapply(p, q)
	if(any(is.nan(y))) warning("NaNs produced")
	return(y)
}
