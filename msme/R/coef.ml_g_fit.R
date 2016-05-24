coef.ml_g_fit <- function(object, ...) {
  object$beta.hat[-length(object$beta.hat)]
}

