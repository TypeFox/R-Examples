powerbasis <- function(x, exponents, nderiv=0) {
#POWERBASIS computes values of monomials, or their derivatives.
#  The powers of X are the NBASIS nonnegative integers in EXPONENTS.
#  The default is 1, meaning X itself.
#  Arguments are as follows:
#  X         ... vector of values at which the polynomials are to
#                evaluated
#  EXPONENTS ... vector of exponents
#  NDERIV    ... order of derivative to be returned.
#  Return is:
#  A matrix with length(X) rows and NBASIS columns containing
#    the values of the monomials or their derivatives

#  last modified 13 December 2002

	x <- as.vector(x)
	n <- length(x)

	nbasis <- length(exponents)

	powermat <- matrix(0,n,nbasis)
	if (nderiv == 0) {
    	for (ibasis in 1:nbasis)
        	powermat[,ibasis] <- x^exponents[ibasis]
	} else {
    	if (any(exponents - nderiv < 0) && any(x == 0)) {
        	stop("A negative exponent is needed and an argument value is 0.")
    	} else {
        	for (ibasis in 1:nbasis) {
            	degree <- exponents[ibasis]
            	if (nderiv <= degree) {
                	fac <- degree
                	for (ideriv in 2:nderiv) {
                    	fac <- fac*(degree-ideriv+1)
                	}
                	powermat[,ibasis] <- fac*x^(degree-nderiv)
            	}
        	}
    	}
	}
	return(powermat)
}
