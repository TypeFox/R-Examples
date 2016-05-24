sqrm <- function(x, ...){
  ## checking if matrix is quadratic
  stopifnot(ncol(x) == nrow(x))
  ## checking if matrix is scalar
  if(ncol(x) == 1){
    ans <- sqrt(x)
    return(ans)
  }
  ## computing square root of matrix
  e <- eigen(x, ...)
  V <- e$vectors
  ans <- V %*% diag(sqrt(e$values)) %*% t(V)
  if(!is.null(dimnames(x))) dimnames(ans) <- dimnames(x)
  return(ans)
}
