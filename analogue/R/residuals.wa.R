`residuals.wa` <- function(object, ...) {
  resi <- object$residuals
  if (!is.null(object$na.action))
    resi <- naresid(object$na.action, resi)
  resi
}
