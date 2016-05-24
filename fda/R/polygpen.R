polygpen <- function(basisobj, Lfdobj=int2Lfd(1))
{

#  Computes the polygonal penalty matrix.
#  Arguments:
#  BASISOBJ ... a basis object of type "polyg"
#  LFDOBJ   ... either the order of derivative or a
#               linear differential operator to be penalized.
#          The highest derivative must be either 0 or 1.
#  Returns the penalty matrix.

#  Last modified 3 January 2008

if (!(inherits(basisobj, "basisfd"))) stop(
    "First argument is not a basis object.")

Lfdobj <- int2Lfd(Lfdobj)

type <- basisobj$type
if (type != "polyg") stop("BASISOBJ not of type polyg")

nderiv <- Lfdobj$nderiv

if (nderiv > 1) stop(
    "Derivative greater than 1 cannot be taken for polygonal basis.")

bwtlist <- Lfdobj$bwtlist
isintLfd <- TRUE
if (nderiv > 0) {
  for (ideriv in 1:nderiv) {
    fdj <- bwtlist[[ideriv]]
    if (!is.null(fdj)) {
      if (any(fdj$coefs != 0)) {
        isintLfd <- FALSE
        break
      }
    }
  }
}

if (isintLfd) {
    args    <- basisobj$params
    n       <- length(args)
    argdiff <- diff(args)
    penaltymat <- diag(rep(1,n))
    if (nderiv == 0) {
        penaltymat[1,1] = argdiff[  1]/3
        penaltymat[n,n] = argdiff[n-1]/3
        indx = 2:(n-1)
        for (i in indx)
            penaltymat[i,i] = (argdiff[i]+argdiff[i-1])/3
        indx = 2:n
        for (i in indx) {
            penaltymat[i,i-1] = argdiff[i-1]/6
            penaltymat[i-1,i] = argdiff[i-1]/6
        }
    } else {
        argdiff = 1/argdiff
        penaltymat[1,1] = argdiff[  1]
        penaltymat[n,n] = argdiff[n-1]
        indx = 2:(n-1)
        for (i in indx)
            penaltymat[i,i] = argdiff[i]+argdiff[i-1]
        indx = 2:n
        for (i in indx) {
            penaltymat[i,i-1] = -argdiff[i-1]
            penaltymat[i-1,i] = -argdiff[i-1]
        }
    }
    penaltymatrix <- penaltymat
} else {
    penaltymatrix <- inprod(basisobj, basisobj, Lfdobj, Lfdobj)
}

return( penaltymatrix )
}
