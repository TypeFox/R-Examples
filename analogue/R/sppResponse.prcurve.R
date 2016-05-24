`sppResponse.prcurve` <- function(x, n = 100, ...) {
  Pred <- function(x, n) {
    rn <- range(x$lambda)
    newx <- seq(from = rn[1], to = rn[2], length = n)
    fitted.values <- predict(x$model, x = newx)
    observed <- x[c("lambda","x")]
    names(fitted.values) <- names(observed) <- c("gradient","response")
    list(observed = observed, fitted.values = fitted.values)
  }
  out <- lapply(x$smooths, Pred, n = n)
  class(out) <- "sppResponse"
  out
}
