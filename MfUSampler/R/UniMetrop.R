MfU.UniMet <- function(x, f, sigma, ...) {
  x.prop <- rnorm(1, x, sigma)
  p <- min(1.0, exp(f(x.prop, ...) - f(x, ...)))
  if (runif(1) < p) {
    ret <- x.prop
    attr(ret, "accept") <- TRUE
  } else {
    ret <- x
    attr(ret, "accept") <- FALSE
  }
  return (ret)
}