# ========================================================================
# %*% - matrix product
# ========================================================================
# DEBUG: how to make this work with outer()?
#

"%*%.rv" <- function(x, y) { ## CHECK: TTRY if is.constant(x) => normal
  return("%**%"(x, y))
}

"%**%" <- function(x, y) { ## CHECK: TTRY if is.constant(x) => normal
  if (! is.rv(x) && ! is.rv(y)) {
    return(.Primitive("%*%")(x, y))
  }
  d <- dim(y)
  if (! is.rv(x) && (is.null(d)) || (length(d) == 2 && d[2] == 1)) {
    n.sims <- .Primitive("max")(rvnsims(x), rvnsims(y), na.rm=FALSE)
    ysim <- sims(as.rv(y), dimensions=TRUE, n.sims=n.sims)
    # Typical case: constant matrix times a rv vector
    AB <- t(.Primitive("%*%")(x, t(ysim)))
    rvsims(AB)
  } else {
    rvmapply(base::crossprod, x=t(as.rv(x)), y=as.rv(y))
  }
}




