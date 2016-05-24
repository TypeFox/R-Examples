
##==============================================================================
## ldp         : Solves Least Distance Programming
##==============================================================================

ldp <- function(G, H, tol=sqrt(.Machine$double.eps),   verbose=TRUE)  {

  ##------------------------------------------------------------------------
  ## 0. Setup problem
  ##------------------------------------------------------------------------
  ## input consistency
  if (! is.matrix(G) & ! is.null(G)) G <- t(as.matrix(G))

  if (is.null(tol)) tol <- sqrt(.Machine$double.eps)

  ## Problem dimension
  Nx     <- ncol(G)   # number of unknowns
  Nin    <- nrow(G)   # number of inequalities
  if (length(H) != Nin)
    stop("cannot solve least distance problem - G and H not compatible")

  IsError <- FALSE
  NW      <- (Nx+1)*(Nin+2) +2*Nin
  storage.mode(G) <- storage.mode(H) <- "double"
  
  sol  <-.Fortran("ldp",G=G,H=H,
         NUnknowns=as.integer(Nx),NConstraints=as.integer(Nin),
         NW=as.integer(NW),X=as.vector(rep(0,Nx)),XNorm=0.,
         W=as.double(rep(0.,NW)),xIndex=as.integer(rep(0,Nin)),
         Mode=as.integer(0),
         verbose=as.logical(verbose),IsError=as.logical(IsError))
  IsError<-sol$IsError 
  
  ## The solution
  X    <- sol$X

  ## Residual of the inequalities
  residual  <- 0
  if (!is.null(G)) {
    ineq     <- G %*% X - H
    residual  <- -sum(ineq[ineq<0])
  }
  ## The solution norm
  solution <- sum ((X)^2)

  xnames <- colnames(G)
  names (X) <- xnames

  return(list(X=X,                            # vector containing the solution of the least distance problem.
              residualNorm=residual,          # scalar, the sum of residuals of violated inequalities
              solutionNorm=solution,          # scalar, the value of the quadratic function at the solution
              IsError=IsError,                # if an error occurred
              type="ldp"))

}
