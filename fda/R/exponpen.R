exponpen <- function(basisobj, Lfdobj=int2Lfd(2))
{

#  Computes the Exponential penalty matrix.
#  Argument:
#  BASISOBJ ... a basis.fd object of type "expon"
#  LFDOBJ   ... either the order of derivative or a
#                linear differential operator to be penalized.
#  Returns the penalty matrix.

#  Last modified 9 February 2007

#  Check BASISOBJ

if (!(inherits(basisobj, "basisfd"))) stop(
    "First argument is not a basis object.")

type <- basisobj$type
if (type != "expon") stop ("Wrong basis type")

#  Check LFDOBJ

Lfdobj <- int2Lfd(Lfdobj)

#  Compute penalty matrix

if (is.integerLfd(Lfdobj)) {
    nderiv  <- Lfdobj$nderiv
    ratevec <- basisobj$params
    nrate   <- length(ratevec)
    penaltymatrix <- matrix(0,nrate,nrate)
    tl <- basisobj$rangeval[1]
    tu <- basisobj$rangeval[2]
    for (irate in 1:nrate) {
      	ratei <- ratevec[irate]
      	for (jrate in 1:irate) {
        	ratej <- ratevec[jrate]
        	ratesum <- ratei + ratej
        	if (ratesum != 0) {
          		penaltymatrix[irate,jrate] <- (ratei*ratej)^nderiv *
              	(exp(ratesum*tu) - exp(ratesum*tl)) / ratesum
        	} else {
          		if (nderiv == 0) penaltymatrix[irate,jrate] <- tu - tl
        	}
        	penaltymatrix[jrate,irate] <- penaltymatrix[irate,jrate]
      	}
	}
} else {
    penaltymatrix <- inprod(basisobj, basisobj, Lfdobj, Lfdobj)
}

penaltymatrix
}
