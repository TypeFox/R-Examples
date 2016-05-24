coefficients.mvdareg <- function (object, ncomp = object$ncomp, conf = .95, ...) {
  if (is.null(object$validation$coefficients)) {
    {
      if (missing(ncomp) || is.null(ncomp)) {
        B <- object$coefficients[, 1:object$ncomp]
      }
      else {
        B <- object$coefficients[, ncomp]
      }
      return(B)
    } 
  } else {
    coefficients.boots(object, ncomp = ncomp, conf = conf)
  }
}