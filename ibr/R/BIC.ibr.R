BIC.ibr <- function(object, ...) {
  r <- object$residuals
  n <- length(r)
  stderr <- sqrt(sum(r^2)/(n-object$finaldf))
  return(log(stderr^2)+log(n)*object$finaldf/n)
}
