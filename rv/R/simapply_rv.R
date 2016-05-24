# ========================================================================
# simapply  -  apply a (numeric) function to the simulations, rowwise, with dimensions
# ========================================================================
# Vectorizes over the simulations of one single rv; 
# for vectorization over a group of rvs, see 'rvmapply', 'rvvapply'.
#

simapply <- function(x, FUN, ...) {
  # Works pretty similarly as rvmapply does
  L <- .sims.as.list(x)
  Args <- .Primitive("c")(list(FUN=FUN, SIMPLIFY=FALSE, USE.NAMES=FALSE), list(L))
  Args$MoreArgs <- list(...)
  S <- do.call(base:::mapply, Args)
  r <- rvsims(S) 
  if (isTRUE(all.equal(dim(r), dim(x)))) {
    dimnames(r) <- dimnames(x)
  }
  return(r)
}
