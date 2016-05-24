project.basis <- function(y, argvals, basisobj, penalize=FALSE, 
                          returnMatrix=FALSE)
{
#  Arguments for this function:
#
#  Y        ... an array containing values of curves
#               If the array is a matrix, rows must correspond to argument
#               values and columns to replications, and it will be assumed
#               that there is only one variable per observation.
#               If Y is a three-dimensional array, the first dimension
#               corresponds to argument values, the second to replications,
#               and the third to variables within replications.
#               If Y is a vector, only one replicate and variable are assumed.
#  ARGVALS  ... A vector of argument values.  This must be of length
#    length(Y) if Y is a vector or dim(Y)[1] otherwise.
#  BASISOBJ ... A basis.fd object
#  PENALIZE ... If TRUE, a penalty term is used to deal with a singular
#               basis matrix.  But this is not normally needed.
#
#  Returns a coefficient vector or array. The first dimension is the number
#     of basis functions and the other dimensions (if any) match
#  the other dimensions of Y.
#
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.

#  Last modified 9 May 2012  by Jim Ramsay

#  Check BASISOBJ

  if (!inherits(basisobj, "basisfd")) stop(
    "BASISOBJ is not a basis object.")

#
#  Calculate the basis and penalty matrices, using the default
#   for the number of derivatives in the penalty.
#
  basismat <- getbasismatrix(argvals, basisobj, 0, returnMatrix)
  Bmat     <- crossprod(basismat)
  if (penalize) {
    penmat <- eval.penalty(basisobj)
#
#  Add a very small multiple of the identity to penmat
#   and find a regularization parameter
#
    penmat <- penmat + 1e-10 * max(penmat) * diag(dim(penmat)[1])
    lambda <- (0.0001 * sum(diag(Bmat)))/sum(diag(penmat))
    Cmat <- Bmat + lambda * penmat
  } else {
    Cmat <- Bmat
  }
#
#  Do the fitting by a simple solution of the
#    equations taking into account smoothing
#
  if (is.array(y) == FALSE) y <- as.array(y)
  if(length(dim(y)) <= 2) {
    Dmat <- crossprod(basismat, y)
    coef <- symsolve(Cmat, Dmat)
  } else {
    nvar <- dim(y)[3]
    coef <- array(0, c(basisobj$nbasis, dim(y)[2], nvar))
    for(ivar in 1:nvar) {
      Dmat <- crossprod(basismat, y[,  , ivar])
      coef[,  , ivar] <- symsolve(Cmat, Dmat)
    }
  }
  coef
}
