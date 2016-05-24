

rep.rv <- function (x, times, ...) {
  if (! is.rv(x)) return(rep(x, times, ...))
  if (missing(times)) {
    a <- rep(unclass(x), ...)
  } else {
    a <- rep(unclass(x), times, ...)
  }
  class(a) <- class(x)
  return(a)
}
