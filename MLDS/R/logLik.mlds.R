`logLik.mlds` <-
function(object, ...) {
#object, obj of class mlds
	if (object$method == "glm")
		val <- logLik(object$obj) else
		{ val <- object$logLik	
		attr(val, "df") <- if(object$method == "optim")
		 	length(object$pscale) - 1 else
		 	length(object$par) + 1
		class(val) <- "logLik"
		}
    val
}

`logLik.mlbs` <- function(object, ...) {
	if (object$method == "glm")
		val <- logLik(object$obj) else
		{val <- object$logLik
		attr(val, "df") <- length(object$par) + 1 		}
	    val
}

