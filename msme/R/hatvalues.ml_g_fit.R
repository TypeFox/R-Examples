hatvalues.ml_g_fit <- function(model, ...) {
  tcrossprod(qr.Q(qr(model$X)))
}
