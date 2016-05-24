residuals.bayescomm <-
function (object, ...) {
  object$call$Y - apply(pnorm(object$trace$z), c(2, 3), mean)
}