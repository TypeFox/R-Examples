fourierpen <- function(basisobj, Lfdobj=int2Lfd(2))
{

  #  Computes the Fourier penalty matrix.
  #  Arguments:
  #  BASISOBJ ... a basis object of type "fourier"
  #  LFDOBJ   ... either the order of derivative or a
  #                linear differential operator to be penalized.
  #  Returns  the penalty matrix.

  #  Note:  The number of basis functions is always odd.  If BASISOBJ
  #  specifies an even number of basis functions, then the number of basis
  #  functions is increased by 1, and this function returns a matrix of
  #  order one larger.

  #  Last modified 9 February 2007

  if (!(inherits(basisobj, "basisfd"))) stop(
		"First argument is not a basis object.")

  nbasis <- basisobj$nbasis
  if (2*(nbasis %/% 2) == nbasis) basisobj$nbasis <- nbasis + 1

  type <- basisobj$type
  if (type != "fourier") stop ("Wrong basis type")

  Lfdobj=int2Lfd(Lfdobj)

  width  <- basisobj$rangeval[2] - basisobj$rangeval[1]
  period <- basisobj$params[1]
  ratio  <- round(width/period)
  nderiv <- Lfdobj$nderiv

  if ((width/period) == ratio && is.integerLfd(Lfdobj)) {

    #  Compute penalty matrix for penalizing integral over one period.

    penaltymatrix <- diag(pendiagfn(basisobj, nderiv))

  } else {

    #  Compute penalty matrix by numerical integration

    penaltymatrix <- inprod(basisobj, basisobj, Lfdobj, Lfdobj)

  }

  return( penaltymatrix )
}

#  ------------------------------------------------------------------

pendiagfn <- function(basisobj, nderiv) {

    nbasis  <- basisobj$nbasis
    period  <- basisobj$params[1]
    rangev  <- basisobj$rangeval
    omega   <- 2*pi/period
    halfper <- period/2
    twonde  <- 2*nderiv
    pendiag <- rep(0,nbasis)
    if (nderiv == 0) pendiag[1] <- period/2.0 else pendiag[1] <- 0
    j   <- seq(2,nbasis-1,2)
    fac <- halfper*(j*omega/2)^twonde
    pendiag[j]   <- fac
    pendiag[j+1] <- fac
    pendiag <- 2*pendiag/period
    return(pendiag)
}
