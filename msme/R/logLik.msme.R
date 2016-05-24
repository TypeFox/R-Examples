logLik.msme <- function(object, ...) {
  val <- object$fit$value
  attr(val, "nall") <- nrow(object$X)
  attr(val, "nobs") <- nrow(object$X)
  attr(val, "df") <- length(object$fit$par)
  class(val) <- "logLik"  
  val
}
