

rvfactor <- function (x, ...) {
  UseMethod("rvfactor")
}

rvfactor.default <- function (x, levels=NULL, ...) {
  f <- as.factor(x)
  a <- sims(as.rv(as.integer(f)))
  rvf <- rvsims(a)
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- if (is.null(levels)) {
    levels(f)
  } else {
    levels
  }
  return(rvf)
}


rvfactor.rv <- function (x, levels=NULL, ...) {
  a <- sims(x)
  f <- as.factor(a)
  levs <- levels(f)
  rvf <- rvsims(array(as.integer(f), dim(a)))
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- if (is.null(levels)) {
    levs
  } else {
    levels
  }
  return(rvf)
}

is.numeric.rvfactor <- function (x) {
  FALSE
}

is.rvfactor <- function (x) {
  UseMethod("is.rvfactor")
}

is.rvfactor.rvfactor <- function (x) {
  TRUE
} 

is.rvfactor.rv <- function (x) {
  all(rvsimapply(x, is.factor))
} 

is.rvfactor.default <- function (x) {
  FALSE
} 

as.rvfactor <- function (x, ...)
{
  if (is.rvfactor(x)) x else rvfactor(x)
} 


as.rv.rvfactor <- function (x, ...) {
  return(x)
  attr(x, "levels") <- NULL
  clx <- class(x)
  clx <- clx[clx!="rvfactor"]
  class(x) <- clx
  return(x)
}

"[.rvfactor" <- function (x, ..., drop = FALSE) {
  y <- NextMethod("[")
  attr(y, "levels") <- attr(x, "levels")
  class(y) <- oldClass(x)
  lev <- levels(x)
  if (drop) {
    exclude <- if (any(is.na(levels(x)))) { NULL } else { NA }
    y <- factor(y, exclude=exclude)
  }
  return(y)
}
