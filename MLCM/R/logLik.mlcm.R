`logLik.mlcm` <-
function(object, ...) {
	if (object$method == "glm")
		val <- logLik(object$obj) else {
		val <- object$logLik
		attr(val, "df") <- length(object$par)
		class(val) <- "logLik"
		}
	val
	}

