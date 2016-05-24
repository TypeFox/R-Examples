qpoilog<-function(p, mu, sig, lower.tail = TRUE, log.p = FALSE){
	if (log.p) p <- exp(p)
	if(!lower.tail) p <- 1 - p
	q <- function(p) suppressWarnings(qfinder("poilog", p, list(mu=mu, sig=sig)))
	y <- sapply(p, q)
	if(any(is.nan(y))) warning("NaNs produced")
	return(y)
}
