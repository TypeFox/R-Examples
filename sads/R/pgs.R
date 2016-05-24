pgs <- function(q, k, S, lower.tail=TRUE, log.p = FALSE) {
	if (any(!is.wholenumber(q))) warning("non integer values in q")
	y <- suppressWarnings(cumsumW("gs", q, list(k=k, S=S), lower.tail, log.p, pad=TRUE))
	if(any(is.nan(y))) warning("NaNs produced")
	return(y)
}
