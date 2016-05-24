hatvalues.msme <- function(model, ...) {
  tcrossprod(qr.Q(qr(model$X)))
}
