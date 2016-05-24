vcov.aodml <- function(object, ...) {
  
	b <- coef(object)
	nbb <- length(b)
	v <- as.matrix(object$varparam[seq(nbb), seq(nbb)])
	nam <- names(b)
	dimnames(v) <- list(nam, nam)
	v

}

vcov.aodql <- function(object, ...) vcov(object$fm, ...)