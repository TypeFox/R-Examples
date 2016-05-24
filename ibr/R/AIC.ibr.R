AIC.ibr <- function(object, ..., k = 2) {
  r <- object$residuals
  n <- length(r)
  stderr <- sqrt(sum(r^2)/(n-object$finaldf))
  return(log(stderr^2)+2*object$finaldf/n)
}

