pls <- function(q, N, alpha, lower.tail=TRUE, log.p=FALSE) {
	if (any(!is.wholenumber(q))) warning("non integer values in q")
	y <- suppressWarnings(cumsumW("ls", q, list(N=N, alpha=alpha), lower.tail, log.p, pad=TRUE))
	if(any(is.nan(y))) warning("NaNs produced")
	return(y)
}
