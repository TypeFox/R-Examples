# CHNOSZ/util.matrix.R
# miscellaneous functions for matrices

invertible.combs <- function(A, nmax=20) {
  # given a matrix A, with number of rows equal to or greater 
  # than the number of columns, return the combinations of row 
  # numbers that constitute invertible square matrices  20111126 jmd
  # start with all combinations
  # nmax: don't consider more than this many rows (to save time for large systems)
  nr <- min(nrow(A), nmax)
  combs <- t(combn(nr, ncol(A)))
  ic <- numeric()
  zero <- sqrt(.Machine$double.eps)
  # loop over the combinations
  for(i in 1:nrow(combs)) {
    # a slow method, actually calculating the inverse
    #tt <- try(solve(A[combs[i,],]),silent=TRUE)
    #if(class(tt)!="try-error") ic <- c(ic,i)
    # it's faster just to test if the determinant is non-zero
    d <- det(A[combs[i,],])
    if(abs(d) > zero) ic <- c(ic, i)
  }
  return(combs[ic,,drop=FALSE])
}

