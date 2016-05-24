loadings.mvdareg <- function (object, ncomp = object$ncomp, conf = .95, ...) {
  if (is.null(object$validation$loadings)) {
    {
      if (missing(ncomp) || is.null(ncomp)) {
        B <- object$loadings[, 1:object$ncomp]
      }
      else {
        B <- object$loadings[, ncomp]
      }
      return(B)
    } 
  } else {
    loadings.boots(object, ncomp = ncomp, conf = conf)
  }
}