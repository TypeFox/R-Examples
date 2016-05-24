# From SamplerCompare, (c) 2010 Madeleine Thompson

# Computes Q such that Q^T Q = R^T R + x x^T.

chud <- function(R, x) {
  p <- nrow(R)
  stopifnot(is.matrix(R) && p==ncol(R))
  stopifnot(is.numeric(x) && length(x)==p)
  L <- .Fortran(dchud, R, p, p, x, 0, 0, 0, 0, 0, numeric(p), numeric(p))
  return(L[[1]])
}

# Computes Q such that Q^T Q = R^T R - x x^T.

chdd <- function(R, x) {
  p <- as.integer(nrow(R))
  z <- as.integer(0)
  R <- as.matrix(R)
  x <- as.numeric(x)
  stopifnot(p==ncol(R) && p==length(x))
  L <- .Fortran(dchdd, R, p, p, x, z, z, z, z, z, numeric(p), numeric(p),
                integer(1))
  info <- L[[12]]
  if (info==-1)
    stop("downdating produced a non-positive-definite matrix")
  return(L[[1]])
}
