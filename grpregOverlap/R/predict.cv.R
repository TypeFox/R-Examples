
## function: predict cv.grpregOverlap
# -------------------------------------------------------------------------------
predict.cv.grpregOverlap <- function(object, X, type=c("link", "response", "class", "coefficients", "vars", "groups", "nvars", "ngroups", "norm"),
                                      latent = FALSE, lambda = object$lambda.min,
                                      which=object$min, ...) {
  type <- match.arg(type)
  predict(object$fit, X=X, type=type, latent=latent, lambda=lambda,
          which=which, ...)
}

coef.cv.grpregOverlap <- function(object, latent = FALSE, lambda = object$lambda.min,
                                   which = object$min, ...) {
  coef(object$fit, lambda=lambda, latent=latent, which=which, ...)
}
# -------------------------------------------------------------------------------
