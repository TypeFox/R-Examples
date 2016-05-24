# Function from package 'fda' (c) 2014

inprod.bspline <- function(fdobj1, fdobj2=fdobj1, nderiv1=0, nderiv2=0) {
#INPROD_BSPLINE  computes matrix of inner products of the derivatives
#  of order DERIV1 and DERIV2 of two functional data objects
#  FD1 and FD2, respectively.
#  These must both have Bspline bases, and these bases must have
#  a common break point sequence.   However, the orders can differ.
#  If only the first argument is present, the inner products of
#  FD1 with itself is taken.  If argument DERIV is not supplied,
#  it is taken to be 0.
#

#  Last modified 26 October 2005

if(!inherits(fdobj1,"fd")) stop("FD1 is not a functional data object.")
if(!inherits(fdobj2,"fd")) stop("FD2 is not a functional data object.")

basis1 <- fdobj1$basis
type1  <- basis1$type
if(type1!="bspline")   stop("FDOBJ1 does not have a B-spline basis.")

range1  <- basis1$rangeval
breaks1 <- c(range1[1], basis1$params, range1[2])
nbasis1 <- basis1$nbasis
norder1 <- nbasis1 - length(breaks1) + 2

basis2 <- fdobj2$basis
type2  <- basis2$type
if(type2!="bspline")    stop("FDOBJ2 does not have a B-spline basis.")

range2  <- basis2$rangeval
breaks2 <- c(range2[1], basis2$params, range2[2])
nbasis2 <- basis2$nbasis
norder2 <- nbasis2 - length(breaks2) + 2

if(any((range1 - range2) != 0)) stop(
	"The argument ranges for FDOBJ1 and FDOBJ2 are not identical.")

#  check that break values are equal and set up common array

if(length(breaks1) != length(breaks2)) stop(
	"The numbers of knots for FDOBJ1 and FDOBJ2 are not identical")

if(any((breaks1 - breaks2) != 0)) stop(
	"The knots for FDOBJ1 and FDOBJ2 are not identical.")
else    breaks <- breaks1

if(length(breaks) < 2) stop(
	"The length of argument BREAKS is less than 2.")

breakdiff <- diff(breaks)
if(min(breakdiff) <= 0) stop(
	"Argument BREAKS is not strictly increasing.")

#  set up the two coefficient matrices

coef1 <- as.matrix(fdobj1$coefs)
coef2 <- as.matrix(fdobj2$coefs)
if(length(dim(coef1)) > 2) stop("FDOBJ1 is not univariate.")
if(length(dim(coef2)) > 2) stop("FDOBJ2 is not univariate.")

nbreaks   <- length(breaks)
ninterval <- nbreaks - 1
nbasis1   <- ninterval + norder1 - 1
nbasis2   <- ninterval + norder2 - 1
if (dim(coef1)[1] != nbasis1 || dim(coef2)[1] != nbasis2)
    stop(paste("Error: coef1 should have length no. breaks1+norder1-2",
           "and coef2 no. breaks2+norder2-2."))

breaks1 <- breaks[1]
breaksn <- breaks[nbreaks]

# The knot sequences are built so that there are no continuity conditions
# at the first and last breaks.  There are k-1 continuity conditions at
# the other breaks.

temp   <- breaks[2:(nbreaks-1)]
knots1 <- c(breaks1*rep(1,norder1),temp,breaksn*rep(1,norder1))
knots2 <- c(breaks1*rep(1,norder2),temp,breaksn*rep(1,norder2))

# Construct  the piecewise polynomial representation of
#    f^(DERIV1) and g^(DERIV2)

nrep1 <- dim(coef1)[2]
polycoef1 <- array(0,c(ninterval,norder1-nderiv1,nrep1))
for (i in 1:nbasis1) {
    #  compute polynomial representation of B(i,norder1,knots1)(x)
    ppBlist <- ppBspline(knots1[i:(i+norder1)])
    Coeff <- ppBlist[[1]]
    index <- ppBlist[[2]]
    # convert the index of the breaks in knots1 to the index in the
    # variable "breaks"
    index <- index + i - norder1
    CoeffD <- ppderiv(Coeff,nderiv1)  # differentiate B(i,norder1,knots1)(x)
    # add the polynomial representation of B(i,norder1,knots1)(x) to f
    if (nrep1 == 1)
       polycoef1[index,,1] <- coef1[i]*CoeffD +
                                        polycoef1[index,,1]
    else {
        for (j in 1:length(index)){
	         temp <- outer(CoeffD[j,],coef1[i,])
	         polycoef1[index[j],,] <-  temp + polycoef1[index[j],,]
        }
    }
}

nrep2 <- dim(coef2)[2]
polycoef2 <- array(0,c(ninterval,norder2-nderiv2,nrep2))
for(i in 1:nbasis2){
    #  compute polynomial representation of B(i,norder2,knots2)(x)
    ppBlist <- ppBspline(knots2[i:(i+norder2)])
    Coeff <- ppBlist[[1]]
    index <- ppBlist[[2]]
    # convert the index of the breaks in knots2 to the index in the
    # variable "breaks"
    index <- index + i - norder2
    CoeffD <- ppderiv(Coeff, nderiv2)  # differentiate B(i,norder2,knots2)(x)
    # add the polynomial representation of B(i,norder2,knots2)(x) to g
    if (nrep2 == 1) polycoef2[index,,1] <- coef2[i]*CoeffD +
                                           polycoef2[index,,1]
    else{
        for(j in 1:length(index)){
	         polycoef2[index[j],,] <- outer(CoeffD[j,],coef2[i,]) +
	                                    polycoef2[index[j],,]
        }
    }
}

# Compute the scalar product between f and g

prodmat <- matrix(0,nrep1,nrep2)
for (j in 1:ninterval){
    # multiply f(i1) and g(i2) piecewise and integrate
    c1 <- as.matrix(polycoef1[j,,])
    c2 <- as.matrix(polycoef2[j,,])
    polyprodmat <- polyprod(c1,c2)
    # compute the coefficients of the anti-derivative
    N <- dim(polyprodmat)[3]
    delta <- breaks[j+1] - breaks[j]
    power <- delta
    prodmati <- matrix(0,nrep1,nrep2)
    for (i in 1:N){
        prodmati <- prodmati + power*polyprodmat[,,N-i+1]/i
        power    <- power*delta
    }
    # add the integral to s
    prodmat <- prodmat + prodmati
}

prodmat

}
