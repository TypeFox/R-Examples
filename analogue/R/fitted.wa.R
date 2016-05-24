`fitted.wa` <- function(object, ...) {
  fits <- object$fitted.values
  if (!is.null(object$na.action))
    fits <- napredict(object$na.action, fits)
  fits
}
