unique.circular <- function (x, ...) {
  z <- unique.default(x, ...)
  circularp(z) <- circularp(x)
  class(z) <- class(x)
  return(z)
}
