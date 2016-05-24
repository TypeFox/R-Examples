# logLik method for loglm objects, to allow use of AIC() and BIC()
# with MASS::loglm, giving comparable results to the use of these
# functions with glm(..., family=poisson) models.

# allow for non-integer frequencies
# allow for zero frequencies, with a zero= argument

logLik.loglm <- function(object, ..., zero=1E-10) {
	fr <- if(!is.null(object$frequencies)) unclass(object$frequencies) else {
				unclass(update(object, keep.frequencies = TRUE)$frequencies)
			}
	df <- prod(dim(fr)) - object$df
	if (any(fr==0)) {
		fr <- as.vector(fr)
		fr[fr==0] <- zero
	}
	structure(sum((log(fr) - 1) * fr - lgamma(fr + 1))  - object$deviance/2,
			df = df, class = "logLik")
}


