predict.cv.hqreg <- function(object, X, lambda = c("lambda.1se","lambda.min"), type = c("response","coefficients","nvars"), ...) {
  type = match.arg(type)
  if (is.character(lambda)) {
    lambda = match.arg(lambda)
    lambda = object[[lambda]]
  } else if(!is.numeric(lambda)) stop("Invalid lambda")
  predict(object$fit, X, lambda = lambda, type = type, ...)
}

coef.cv.hqreg <- function(object, lambda = c("lambda.1se","lambda.min"), ...) {
  type = match.arg(type)
  if (is.character(lambda)) {
    lambda = match.arg(lambda)
    lambda = object[[lambda]]
  } else if(!is.numeric(lambda)) stop("Invalid lambda")
  coef(object$fit, lambda = lambda, type = type, ...)
}