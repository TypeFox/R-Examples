prbs <- function(q, N, S, lower.tail=TRUE, log.p = FALSE){
	if (any(!is.wholenumber(q))) warning("non integer values in q")
	y <- suppressWarnings(cumsumW("rbs", q, list(N=N, S=S), lower.tail, log.p, pad=TRUE))
	if(any(is.nan(y))) warning("NaNs produced")
	return(y)
}
