logLik.km <- function(object, ...) {

	logLik <-  -0.5*(object@n*log(2*pi) + 2*sum(log(diag(object@T))) + t(object@z)%*%object@z)     
	
	if (object@method == "PMLE") {
		fun <- match.fun(object@penalty.fun)
		param <- object@covariance@range.val
      	penalty <- -object@n * sum(fun(1/param, object@penalty.value))
      logLik <- logLik + penalty
    }
	
	return(logLik)
}

