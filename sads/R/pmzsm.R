pmzsm <- function(q, J, theta, lower.tail=TRUE, log.p=FALSE){
	if (any(!is.wholenumber(q))) warning("non integer values in q")
	y <- suppressWarnings(cumsumW("mzsm", q, list(J=J, theta=theta), lower.tail, log.p, pad=TRUE))
	if(any(is.nan(y))) warning("NaNs produced")
	return(y)
}
