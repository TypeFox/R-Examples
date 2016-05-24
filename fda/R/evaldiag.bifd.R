evaldiag.bifd <- function(evalarg, bifdobj, sLfd=int2Lfd(0), tLfd=int2Lfd(0), 
                            returnMatrix=FALSE)
{
#  EVALDIAG_BIFD  evaluates a bi-functional data object BIFD
#  with both argument values in array EVALARG.
#  SLfd and TLfd are either integers giving the order of derivative,
#  or linear differential operators to be applied before evaluation.
#  Their defaults are 0, meaning that the function itself is evaluated.

#  last modified 2008(?) replacing Matlab subscripts with R style
# in lines 70, 77, 78;  previously modified 26 October 2005
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.

#  Last modified 8 May 2012 by Jim Ramsay

#  exchange order if BIFD is the first argument
if (inherits(evalarg, "bifd")) {
    temp    <- bifdobj
    bifdobj <- evalarg
    evalarg <- temp
}

#  check EVALARG

evalarg <- as.vector(evalarg)

if (!inherits(bifdobj, "bifd")) stop(
    "Argument BIFD is not a bivariate functional data object.")

n <- length(evalarg)

#  extract the two bases

sbasisobj <- bifdobj$sbasis
tbasisobj <- bifdobj$tbasis
snbasis   <- sbasisobj$nbasis
tnbasis   <- tbasisobj$nbasis
ranges    <- sbasisobj$rangeval
ranget    <- tbasisobj$rangeval

#  check that the bases have the same range

if (any(ranges != ranget)) stop(
    "The ranges are not identical.")

#  check the differential operators

sLfd <- int2Lfd(sLfd)
tLfd <- int2Lfd(tLfd)

#  compute the basis matrix for SBASISOBJ

snderiv   <- sLfd$nderiv
sbasismat <- eval.basis(evalarg, sbasisobj, sLfd, returnMatrix)

#  compute the basis matrix for tBASISOBJ

tnderiv   <- tLfd$nderiv
tbasismat <- eval.basis(evalarg, tbasisobj, tLfd, returnMatrix)

#  Extract the coefficient matrix

coef  <- bifdobj$coefs
coefd <- dim(coef)
ndim  <- length(coefd)

if        (ndim == 2) {
        evalarray <- diag(sbasismat %*% coef %*% t(tbasismat))
} else if (ndim == 3) {
        ncurves   <- coefd[3]
        evalarray <- matrix(0,n,ncurves)
        for (i in 1:ncurves)
            evalarray[,i] <- diag(sbasismat %*% coef[,,i] %*% t(tbasismat))
} else if (ndim == 4) {
        ncurves  <- coefd[3]
        nvar     <- coefd[4]
        evalarray <- array(0,c(n,ncurves,nvar))
        for (i in 1:ncurves) {
            for (j in 1:nvar) {
                evalarray[,i,j] <-
                    diag(sbasismat %*% coef[,,i,j] %*% t(tbasismat))
            }
        }
} else {
       stop("The coefficient array has improper dimension.")
}

return(evalarray)
}
