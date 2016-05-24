weights.mvdareg <- function (object, ncomp = object$ncomp, conf = .95, ...) {
  if (is.null(object$validation$weights)) {
    {
      if (missing(ncomp) || is.null(ncomp)) {
        B <- object$weights[, 1:object$ncomp]
      }
      else {
        B <- object$weights[, ncomp]
      }
      return(B)
    } 
  } else {
    weight.boots(object, ncomp = ncomp, conf = conf)
  }
}