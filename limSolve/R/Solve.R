
## =============================================================================
## Solve       : Generalised inverse solution of Ax=B
## =============================================================================

Solve <- function(A, B=diag(nrow=nrow(A)),
   tol=sqrt(.Machine$double.eps))   {

  M <-ginv(A,tol)
  if (is.null(M))
    return(NULL)
  B <- matrix(data=B, nrow=nrow(A))
  X <- M %*% B
  if (ncol(B) == 1)  {
    xnames <- colnames(A)
    X <- as.vector(X)
    names (X) <- xnames
  }
  return(X)
}

