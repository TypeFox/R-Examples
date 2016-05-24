

rvcut <- function (x, ...) {
  UseMethod("rvcut")
}

rvcut.default <- function (x, ...) {
  f <- cut(x, ...)
  levs <- levels(f)
  rvf <- rvsims(as.integer(f))
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- levs
  return(rvf)
}

rvcut.rv <- function (x, ...) {
  a <- sims(x)
  f <- cut(a, ...)
  levs <- levels(f)
  rvf <- rvsims(array(as.integer(f), dim(a)))
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- levs
  return(rvf)
}


