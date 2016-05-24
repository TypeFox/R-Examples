powerpen <- function(basisobj, Lfdobj=int2Lfd(2)) {
#  POWERPEN  Computes the power basis penalty matrix.
#  Arguments:
#  BASISFD  ... a power basis object
#  Lfd      ... either the order of derivative or a
#               linear differential operator to be penalized.
#  Returns a list the first element of which is the basis matrix
#   and the second element of which is the diagonal of the penalty matrix.

#  Last modified:  2008.10.17 by Spencer Graves
#  Previously modified:  3 January 2008 by Jim Ramsay

#  check BASISOBJ

if (!inherits(basisobj, "basisfd")) stop(
    	"First argument is not a basis object.")

if (basisobj$type != "power") stop("BASISOBJ not of type POWER.")

rangeval  <- basisobj$rangeval
exponents <- basisobj$params

#  check LFDOBJ

Lfdobj <- int2Lfd(Lfdobj)

if (is.integerLfd(Lfdobj)) {

    #  case where LFDOBJ is integer

    nderiv <- Lfdobj$nderiv

    if (any(exponents - nderiv < 0) && rangeval[1] <= 0)
        	stop(paste("A negative exponent is needed and",
                    "an argument value can be nonpositive."))
    nbasis     <- basisobj$nbasis
    penaltymat <- matrix(0,nbasis,nbasis)
    xrange     <- basisobj$rangeval
    for (ibasis in 1:nbasis) {
        ideg <- exponents[ibasis]
        if (nderiv == 0) {
            ifac <- 1
        } else {
            ifac <- ideg
            if(nderiv>1) for (k in 2:nderiv) {
                if (ideg == k-1) ifac = 0
                else             ifac <- ifac*(ideg - k + 1)
            }
        }
        if (ibasis > 1) {
            for (jbasis in 1:(ibasis-1)) {
                jdeg <- exponents[jbasis]
        	    if (nderiv == 0) {
                    jfac <- 1
                } else {
                    jfac <- jdeg
                    if(nderiv>1) for (k in 2:nderiv) {
                        if (jdeg == k-1) jfac = 0
                        else             jfac <- jfac*(jdeg - k + 1)
                    }
                }
                if (ifac*jfac == 0) {
                    penaltymat[ibasis,jbasis] = 0
                } else {
                    penaltymat[ibasis,jbasis] <- ifac*jfac*
	                      (xrange[2]^(ideg+jdeg-2*nderiv+1) -
	                       xrange[1]^(ideg+jdeg-2*nderiv+1))/
                            (ideg + jdeg - 2*nderiv + 1)
                }
	          penaltymat[jbasis,ibasis] <- penaltymat[ibasis,jbasis]
            }
        }
        if (ifac == 0)
            penaltymat[ibasis,ibasis] <- 0
        if (2*ideg - 2*nderiv + 1 == 0)
            penaltymat[ibasis,ibasis] <- ifac^2*
                    (log(xrange[2]) - log(xrange[1]))
        if (ifac*(2*ideg - 2*nderiv + 1) != 0)
                penaltymat[ibasis,ibasis] <- ifac^2*
	              (xrange[2]^(2*ideg-2*nderiv+1) -
	               xrange[1]^(2*ideg-2*nderiv+1))/
                    (2*ideg - 2*nderiv + 1)
    }
} else {
    	penaltymat <- inprod(basisobj, basisobj, Lfdobj, Lfdobj)
}

penaltymat

}

