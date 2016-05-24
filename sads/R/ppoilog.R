ppoilog <- function(q, mu, sig, lower.tail=TRUE, log.p=FALSE){
	if (any(!is.wholenumber(q))) warning("non integer values in q")
	y <- suppressWarnings(cumsumW("poilog", q, list(mu=mu, sig=sig), lower.tail, log.p, pad=FALSE))
	if(any(is.nan(y))) warning("NaNs produced")
	return(y)
}
