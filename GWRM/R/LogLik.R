#' @export
logLik.gw <- function (object, ...){
  val <- object$loglik
  attr(val, "nobs") <- object$nobs
  attr(val, "df") <- length(object$coefficients)
  class(val) <- "logLik"
  val
}