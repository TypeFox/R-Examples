lambda2df <- function (argvals, basisobj, wtvec=rep(1,n), Lfdobj=NULL, 
                       lambda=0, returnMatrix=FALSE)
{
  #  Computes the the degrees of freedom associated with a regularized
  #    basis smooth by calculating the trace of the smoothing matrix.

  #  Arguments for this function:
  #
  #  ARGVALS  ... A set of argument values.
  #  BASISOBJ ... A basis.fd object created by function create.basis.fd.
  #  WTVEC    ... A vector of N weights, set to one by default, that can
  #               be used to differentially weight observations in the
  #               smoothing phase
  #  LFDOBJ   ... The order of derivative or a linear differential
  #               operator to be penalized in the smoothing phase.
  #               By default Lfdobj is set in function GETBASISPENALTY
  #  LAMBDA   ... The smoothing parameter determining the weight to be
  #               placed on the size of the derivative in smoothing.  This
  #               is 0 by default.
  #  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
  #               from a call to function BsplineS.  See this function for
  #               enabling this option.
  #  Returns:
  #  DF    ...  a degrees of freedom measure

  #  Last modified:  7 May 2012

  n        <- length(argvals)
  nbasis   <- basisobj$nbasis
  if (lambda == 0) {
    df <- nbasis
    return( df )
  }
  if (length(wtvec) != n) stop("WTVEC of wrong length")
  if (min(wtvec) <= 0)    stop("All values of WTVEC must be positive.")
  basismat <- getbasismatrix(argvals, basisobj, 0, returnMatrix)
  basisw   <- basismat*outer(wtvec,rep(1,nbasis))
  Bmat     <- crossprod(basisw,basismat)
  penmat   <- getbasispenalty(basisobj, Lfdobj)
  Bnorm    <- sqrt(sum(Bmat^2))
  pennorm  <- sqrt(sum(penmat^2))
  condno   <- pennorm/Bnorm
  if (lambda*condno > 1e12) {
    lambda <- 1e12/condno
    warning(paste("lambda reduced to",lambda,"to prevent overflow"))
  }
  Cmat     <- Bmat + lambda*penmat
  Cmat     <- (Cmat + t(Cmat))/2
  if (is.diag(Cmat)) {
      Cmatinv <- diag(1/diag(Cmat))
  } else {
      Lmat    <- chol(Cmat)
      Lmatinv <- solve(Lmat)
      Cmatinv <- crossprod(t(Lmatinv))
  }
  hatmat <- Cmatinv %*% Bmat
  df <- sum(diag(hatmat))
  return( df )
}
