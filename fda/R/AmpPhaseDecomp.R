AmpPhaseDecomp <- function(xfd, yfd, hfd, rng=xrng, returnMatrix=FALSE)
{
#  Computes the amplitude-phase decomposition for a registration.

#  Arguments:
#  XFD  ...  FD object for unregistered functions
#  YFD  ...  FD object for registered functions
#  HFD  ...  FD object for warping functions

#  Returns:
#  MS.amp ... mean square for amplitude variation 
#  MS.pha ... mean square for amplitude variation 
#  RSQR   ... squared correlation measure of prop. phase variation 
#  C      ... constant C
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.

#  Last modified 8 May 2012 by Jim Ramsay

xbasis  <- xfd$basis
nxbasis <- xbasis$nbasis
nfine   <- max(201,10*nxbasis)
xrng    <- xbasis$rangeval

if (rng[1] < xrng[1] || rng[2] > xrng[2]) stop(
    "RNG is not within the range of the other arguments.")
if (rng[1] >= rng[2]) stop("Elements of rng are not increasing.")
tfine   <- seq(rng[1],rng[2],len=nfine)
delta   <- tfine[2] - tfine[1]

Dhfine  <- eval.fd(tfine, hfd, 1, returnMatrix)
xfine   <- eval.fd(tfine, xfd, 0, returnMatrix)
yfine   <- eval.fd(tfine, yfd, 0, returnMatrix)
mufine  <- apply(xfine, 1, mean)
etafine <- apply(yfine, 1, mean)

N       <- dim(xfine)[2]
rfine   <- yfine - outer(etafine,rep(1,N))

intetasqr <- delta*trapz(etafine^2)
intmusqr  <- delta*trapz(mufine^2)

covDhSy <- rep(0,nfine)
for (i in 1:nfine) {
    Dhi        <- Dhfine[i,]
    Syi        <- yfine[i,]^2
    covDhSy[i] <- cov(Dhi, Syi)
}
intcovDhSy <- delta*trapz(covDhSy)

intysqr <- rep(0,N)
intrsqr <- rep(0,N)
for (i in 1:N) {
    intysqr[i] <- delta*trapz(yfine[,i]^2)
    intrsqr[i] <- delta*trapz(rfine[,i]^2)
}

C      <- 1 + intcovDhSy/mean(intysqr)
MS.amp <- C*mean(intrsqr)
MS.pha <- C*intetasqr - intmusqr
RSQR   <- MS.pha/(MS.amp+MS.pha)

return(list("MS.amp" = MS.amp, "MS.pha" = MS.pha, "RSQR" = RSQR, "C" = C)) 

}

trapz = function(x) {
n = length(x)
intx = sum(x) - 0.5*(x[1]+x[n])
return(intx)
}

