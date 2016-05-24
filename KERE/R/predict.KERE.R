predict.KERE <- function(object, kern, x, newx,...) {
  nfit <- kernelMult(kern, newx, x, object$alpha[-1,])
  nfit <- sweep(nfit, MARGIN = 2, object$alpha[1,], "+")
  nfit
}
