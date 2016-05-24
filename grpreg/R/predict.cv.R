predict.cv.grpreg <- function(object, X, lambda=object$lambda.min, which=object$min, type=c("link", "response", "class", "coefficients", "vars", "groups", "", "nvars", "ngroups", "norm"), ...) {
  type <- match.arg(type)
  predict.grpreg(object$fit, X=X, lambda=lambda, which=which, type=type, ...)
}
coef.cv.grpreg <- function(object, lambda=object$lambda.min, which=object$min, ...) {
  coef.grpreg(object$fit, lambda=lambda, which=which, ...)
}
