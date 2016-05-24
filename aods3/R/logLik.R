logLik.aodml <- function(object, ...) {
  
  z <- object$logL
  attr(z, "df") <- object$df.model
  attr(z, "nobs") <- object$df.model + object$df.residual
	class(z) <- "logLik"
	z

}

logLik.aodql <- function(object, ...) logLik(object$fm, ...)


