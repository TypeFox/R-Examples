monomialpen <- function(basisobj, Lfdobj=int2Lfd(2),
                        rng=basisobj$rangeval)
{
#  MONOMIALPEN  Computes a monomial basis penalty matrix.
#
#  Warning:  This version is incomplete in the sense that
#    it only works with LFDOBJ = D^m
#
#  Arguments:
#  BASISOBJ ... a monomial basis object
#  LFDOBJ   ... either the order of derivative or a
#               linear differential operator to be penalized.
#  RNG      ... a range over which the penalty matrix is computed.
#  Returns a list the first element of which is the basis matrix
#   and the second element of which is the diagonal of the penalty matrix.

#  Last modified:  July 6, 2012 by Spencer Graves
#    to support rng of class Date and POSIXct
#  Previously modified:  26 October 2005

#  check BASISOBJ

  if (!inherits(basisobj, "basisfd")) stop(
		"First argument is not a basis.fd object.")

  type <- basisobj$type
  if (!(type == "monom")) stop("BASISOBJ not of type monom")

#  check LFDOBJ

  Lfdobj <- int2Lfd(Lfdobj)

#  check whether LFDOBJ is of form D^m

  nderiv <- Lfdobj$nderiv

#  get basis information

  nbasis    <- basisobj$nbasis
  exponents <- basisobj$params

#  check exponents

  for (ibasis in 1:nbasis) {
	ideg <- exponents[ibasis]
	if (ideg-floor(ideg) != 0) stop(
		"An exponent is not an integer.")
    }

#  check rng
  Rng <- rng
  if(!is.numeric(Rng)){
    op <- options(warn=-1)
    rng <- as.numeric(Rng)
    options(op)
    nNAr <- sum(is.na(rng))
    if(nNAr>0)
      stop('as.numeric(rng) contains ', nNAr,
           ' NA', c('', 's')[1+(nNAr>1)],
           ';  class(rng) = ', class(Rng))
  }

#  compute the penalty matrix

  penaltymat <- matrix(0,nbasis,nbasis)
  for (ibasis in 1:nbasis) {
    ideg <- exponents[ibasis]
    ifac <- 1
    if (nderiv > 0) {
      ifac <- ideg
      if (nderiv > 1) for (k in 2:nderiv) ifac <- ifac*(ideg - k + 1)
    }
    for (jbasis in 1:ibasis) {
      jdeg <- exponents[jbasis]
      jfac <- 1
      if (nderiv > 0) {
        jfac <- jdeg
        if (nderiv > 1) for (k in 2:nderiv) jfac <- jfac*(jdeg - k + 1)
      }
      if ((ideg >= nderiv) & (jdeg >= nderiv)) {
        ipow <- ideg+jdeg-2*nderiv+1
        penaltymat[ibasis,jbasis] <-
            (rng[2]^ipow - rng[1]^ipow)*ifac*jfac/ipow
        penaltymat[jbasis,ibasis] <- penaltymat[ibasis,jbasis]
      }
    }
  }
penaltymat
}
