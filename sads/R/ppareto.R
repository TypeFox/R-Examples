ppareto <- function(q, shape, scale = min(q), lower.tail = TRUE, log.p = FALSE){
	shape[ !is.finite(shape) | shape <= 0] <- NaN
	scale[ !is.finite(scale) | scale <= 0] <- NaN
	y <- 1 - (scale/q)^shape
	y [ q < scale ] <- 0
	if (!lower.tail) y <- 1 - y
	if (log.p) y <- log(y)
	if (any(is.nan(y))) warning ("NaNs produced")
	return(y)
}
