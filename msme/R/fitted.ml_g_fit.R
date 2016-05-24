fitted.ml_g_fit <- function(object, ...) {
   as.numeric(object$X %*% coef(object))
}
