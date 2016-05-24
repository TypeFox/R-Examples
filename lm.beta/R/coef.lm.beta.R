coef.lm.beta <- function(object, standardized=TRUE, ...) {
  if(standardized) {
    res <- object$standardized.coefficients
  } else {
    res <- object$coefficients
  }
  res
}