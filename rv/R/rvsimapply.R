# ========================================================================
# rvsimapply  -  apply a function to the simulations, componentwise
# ========================================================================

rvsimapply <- function(x, FUN, ...)
{
  dx <- dim(x)
  n <- length(x)
  if (n==0) {
    return(NULL)
  }
  mv <- lapply(unclass(x), FUN, ...)
  lmv <- sapply(mv, length)
  if (all(lmv==1)) {
    m <- unlist(mv, use.names=TRUE)
    dim(m) <- dx
    dimnames(m) <- dimnames(x)
    return(m)
  } else if (all(lmv==rvnsims(x))) {
    # simulation-wise function was applied - return an object of same type
    attributes(mv) <- attributes(x)
    return(mv)
  } else if (all(lmv==lmv[1])) {
    m <- unlist(mv)
    m <- matrix(m, nrow=lmv[1], ncol=n)
    dimnames(m) <- list(names(mv[[1]]), names(x))
    return(m)
  } else {
    names(mv) <- names(x)
    return(mv)
  }
}

