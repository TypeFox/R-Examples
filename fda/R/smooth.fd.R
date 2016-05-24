smooth.fd <- function(fdobj, fdParobj){
#SMOOTH_FD smooths a functional data object.
#
#  Arguments for this function:
#
#  FDOBJ    ... A functional data object.
#  FDPAROBJ ... A functional parameter object.
#
#  Returns a functional data object containing a smoothed version
#    of the input functional data object
#
#  Last modified: 20081003;  previously modified 26 October 2005
#

#  check fdParobj

if (!inherits(fdParobj,"fdPar")) stop(
		"FDPAROBJ is not a functional parameter object.")

#  check LFD

Lfdobj <- fdParobj$Lfd
Lfdobj <- int2Lfd(Lfdobj)
nderiv <- Lfdobj$nderiv

#  set up FDOBJ

newfdobj <- fdParobj$fd

#  set up basis

basisobj <- newfdobj$basis

#
#  Main smoothing step
#

coef  <- fdobj$coefs
coefd <- dim(coef)
ndim  <- length(coefd)
if (ndim == 3)  nvar <- coefd[3] else nvar <- 1

Bmat  <- inprod(basisobj, basisobj)

#
#  set up coefficient matrix for normal equations
#

lambda <- fdParobj$lambda
penmat <- eval.penalty(basisobj, Lfdobj)

penmat <- eval.penalty(basisobj, Lfdobj)
Cmat   <- Bmat + lambda * penmat

#
#  solve normal equations for each observation
#
if (ndim < 3){
     Dmat <- inprod(basisobj, fdobj)
     coef <- solve(Cmat, Dmat)
} else {
	coef <- array(0,coefd)
    for(ivar in 1:nvar){
        Dmat <- inprod(basisobj, fdobj[,ivar])
        coef[,,ivar] <- solve(Cmat, Dmat)
    }
}

#  set up the smoothed functional data object

fdnames <- fdobj$fdnames
smthfd  <- fd(coef, basisobj, fdnames)

return(smthfd)
}
