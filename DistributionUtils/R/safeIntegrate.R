safeIntegrate <- function(f, lower, upper, subdivisions = 100,
                          rel.tol = .Machine$double.eps^0.25,
                          abs.tol = rel.tol, stop.on.error = TRUE,
                          keep.xy = FALSE, aux = NULL, ...) {

  f <- match.fun(f)
  ff <- function(x) f(x, ...)
  equal <- all.equal(lower, upper)

  if (is.character(equal))
    equal <- FALSE

  if (equal) {
    if (!is.finite(lower)) {
      res <- NULL
      res$value <- 0
      res$abs.error <- 0
      res$subdivisions <- NA
      res$message <- "OK"
      res$call <- match.call()
      class(res) <- "integrate"
    } else {
      res <- NULL
      res$value <- 0
      res$abs.error <- mean(c(ff(upper),ff(lower)))*(upper - lower)
      res$subdivisions <- NA
      res$message <- "OK"
      res$call <- match.call()
      class(res) <- "integrate"
    }
  } else {
    res <- integrate(ff, lower, upper, ...)
  }

  res
}
