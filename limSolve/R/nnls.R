
##==============================================================================
## nnls        : Solves nonnegative least squares problem
##==============================================================================

nnls <- function(A, B, tol=sqrt(.Machine$double.eps), verbose=TRUE) {

##------------------------------------------------------------------------
## 0. Setup problem
##------------------------------------------------------------------------

  ## input consistency
  if (! is.matrix(A) & ! is.null(A))
    A <- t(as.matrix(A))

  if (is.null(tol))
    tol <- sqrt(.Machine$double.eps)

  ## Problem dimension
  Nx     <- ncol(A)   # number of unknowns
  Neq    <- nrow(A)   # number of inequalities
  if (length(B) != Neq)
    stop("cannot solve nnls problem - A and B not compatible")

  sol  <-.Fortran("xnnls", A = as.double(A), MDA = as.integer(Neq),
                  M = as.integer(Neq),
                  N = as.integer(Nx), B = as.double(B), 
                  X = as.vector(rep(0,Nx)), RNorm = 0.,
                  W = as.double(rep(0.,Nx)),ZZ = as.double(rep(0.,Neq)),
                  Index = as.integer(rep(0,Nx)), Mode = as.integer(0))
  IsError <- FALSE
  Mode <- sol$Mode

  if (Mode != 1)  {
    IsError <- TRUE
    if (verbose) {
      if (Mode ==2)
        warning("nnls: The dimensions of the problem are bad")
      if (Mode ==3)
        warning("nnls: iteration count exceeded: more than 3*N iterations")
    }
  }
  
  ## The solution
  X    <- sol$X

  ## Residual of the inequalities
  residual  <- -sum(X[X<0])

  ## The solution norm
  solution <- sum(abs(A %*% X - B))

  xnames <- colnames(A)
  names (X) <- xnames

  return(list(X=X,
              residualNorm=residual,   # sum of violated inequalities
              solutionNorm=solution,   # the quadratic function
              IsError=IsError,
              type="nnls"))

}
