pvolkov <- function(q, theta, m , J, lower.tail=TRUE, log.p=FALSE){
	if (any(!is.wholenumber(q))) warning("non integer values in q")
	y <- suppressWarnings(cumsumW("volkov", q, list(theta=theta, m=m, J=J), lower.tail, log.p, pad=TRUE))
	if(any(is.nan(y))) warning("NaNs produced")
	return(y)
}
