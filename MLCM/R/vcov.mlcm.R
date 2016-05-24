`vcov.mlcm` <- function(object, ...)
	if (object$method == "glm")
		vcov(object$obj, ...) else
		solve(object$hess)