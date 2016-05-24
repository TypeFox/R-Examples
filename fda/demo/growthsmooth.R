#  The growh data are smoothed using two techniques:
#  1.  a standard non-monotone smoothing using function smooth.basis
#      with a penalty on the 4th derivative in order to estimate a
#      smooth acceleration function.
#  2.  a monotone smoothing using function smooth.monotone

#  Last modified 9 March 2007

#  load the data

load("growthdata")

age  <- growthdata$age
hgtm <- growthdata$hgtm
hgtf <- growthdata$hgtf
nage <- length(age)

#  --------------------------------------------------------------------
#                   Smooth the data non-monotonically  
#  --------------------------------------------------------------------

#  This smooth uses the usual smoothing methods to smooth the data,
#  but is not guaranteed to produce a monotone fit.  This may not
#  matter much for the estimate of the height function, but it can
#  have much more serious consequences for the velocity and
#  accelerations.  See the monotone smoothing method below for a
#  better solution, but one with a much heavier calculation overhead.

#  A B-spline basis with knots at age values and order 6 is used

rng      <- c(1,18)
knots    <- age
norder   <- 6
nbasis   <- length(knots) + norder - 2
hgtbasis <- create.bspline.basis(rng, nbasis, norder, knots)
agefine  <- seq(1,18,length=101)

#  --- Smooth these objects, penalizing the 4th derivative  --
#  This gives a smoother estimate of the acceleration functions

Lfdobj    <- 4
lambda    <- 1e-2
growfdPar <- fdPar(hgtbasis, Lfdobj, lambda)

hgtmfd <- smooth.basis(age, hgtm, growfdPar)$fd
hgtffd <- smooth.basis(age, hgtf, growfdPar)$fd

#  plot data and smooth, residuals, velocity, and acceleration

#  Males:

hgtmfitN <- eval.fd(age,     hgtmfd)
hgtmhatN <- eval.fd(agefine, hgtmfd)
velmhatN <- eval.fd(agefine, hgtmfd, 1)
accmhatN <- eval.fd(agefine, hgtmfd, 2)

par(mfrow=c(2,2),pty="s",ask=T)
children <- 1:ncasem
for (i in children) {
    plot(age, hgtm[,i], ylim=c(60,200),
         xlab="", ylab="cm", main=paste("Height for male",i))
    lines(agefine, hgtmhatN[,i], col=4)
    resi <- hgtm[,i] - hgtmfitN[,i]
    ind  <- resi >= -.7 & resi <= .7
    plot(age[ind], resi[ind], type="b", ylim=c(-.7,.7), col=4,
         xlab="", ylab="cm", main="Residuals")
    abline(h=0, lty=2)
    ind <- velmhatN[,i] >= 0 & velmhatN[,i] <= 20
    plot(agefine[ind], velmhatN[ind,i], type="l", ylim=c(0,20), col=4,
         xlab="Years", ylab="cm/yr", main="Velocity")
    abline(h=0, lty=2)
    ind <- accmhatN[,i] >= -6 & accmhatN[,i] <= 6
    plot(agefine[ind], accmhatN[ind,i], type="l", ylim=c(-6,6), col=4,
         xlab="Years", ylab="cm/yr/yr", main="Acceleration")
    abline(h=0, lty=2)
}    

# Females:

hgtffitN <- eval.fd(age,     hgtffd)
hgtfhatN <- eval.fd(agefine, hgtffd)
velfhatN <- eval.fd(agefine, hgtffd, 1)
accfhatN <- eval.fd(agefine, hgtffd, 2)

par(mfrow=c(2,2),pty="s",ask=T)
children <- 1:ncasef
for (i in children) {
    plot(age, hgtf[,i], ylim=c(60,200),
         xlab="", ylab="cm", main=paste("Height for female",i))
    lines(agefine, hgtfhatN[,i], col=4)
    resi <- hgtf[,i] - hgtffitN[,i]
    ind  <- resi >= -.7 & resi <= .7
    plot(age[ind], resi[ind], type="b", ylim=c(-.7,.7), col=4,
         xlab="", ylab="cm", main="Residuals")
    abline(h=0, lty=2)
    ind <- velfhatN[,i] >= 0 & velfhatN[,i] <= 20
    plot(agefine[ind], velfhatN[ind,i], type="l", ylim=c(0,20), col=4,
         xlab="Years", ylab="cm/yr", main="Velocity")
    abline(h=0, lty=2)
    ind <- accfhatN[,i] >= -6 & accfhatN[,i] <= 6
    plot(agefine[ind], accfhatN[ind,i], type="l", ylim=c(-6,6), col=4,
         xlab="Years", ylab="cm/yr/yr", main="Acceleration")
    abline(h=0, lty=2)
}    

#  -------------------------------------------------------------------
#         save the results of the non-monotone smooths
#  -------------------------------------------------------------------

hgtmfdPar <- fdPar(hgtmfd, Lfdobj, lambda)
hgtffdPar <- fdPar(hgtffd, Lfdobj, lambda)

growthfd <- list(hgtmfdPar  = hgtmfdPar,  hgtffdPar  = hgtffdPar)

save(growthfd, file="growthfd")

#  --------------------------------------------------------------------
#                  Compute monotone smooths of the data  
#  --------------------------------------------------------------------

#  These analyses use a function written entirely in S-PLUS called
#  smooth.monotone that fits the data with a function of the form
#                   f(x) = b_0 + b_1 D^{-1} exp W(x)
#     where  W  is a function defined over the same range as X,
#                 W + ln b_1 = log Df and w = D W = D^2f/Df.
#  The constant term b_0 in turn can be a linear combinations of covariates:
#                         b_0 = zmat * c.
#  The fitting criterion is penalized mean squared error:
#    PENSSE(lambda) = \sum [y_i - f(x_i)]^2 +
#                     \lambda * \int [L W(x)]^2 dx
#  where L is a linear differential operator defined in argument Lfdobj.
#  The function W(x) is expanded by the basis in functional data object
#  Because the fit must be calculated iteratively, and because S-PLUS
#  is so slow with loopy calculations, these fits are VERY slow.  But
#  they are best quality fits that I and my colleagues, notably
#  R. D. Bock, have been able to achieve to date.
#  The Matlab version of this function is much faster.

#  ------  First set up a basis for monotone smooth   --------

#  We use b-spline basis functions of order 6 
#  Knots are positioned at the ages of observation.

rng     <- c(1,18)
norder  <- 6
nbasis  <- nage + norder - 2
wbasis  <- create.bspline.basis(rng, nbasis, norder, age)
agefine <- seq(1,18,length=101)

wgt <- rep(1,nage)

#  starting values for coefficient

cvec0 <- matrix(0,nbasis,1)
Wfd0  <- fd(cvec0, wbasis)

Lfdobj    <- 3          #  penalize curvature of acceleration
lambda    <- 10^(-1)  #  smoothing parameter
growfdPar <- fdPar(Wfd0, Lfdobj, lambda)

#  ---------------------  Now smooth the data  --------------------

# Males:

cvecm <- matrix(0, nbasis, ncasem)
betam <- matrix(0, 2,      ncasem)
RMSEm <- matrix(0, 1,      ncasem)

#  setting the output to unbuffered mode in Misc menu item might
#  be appreciated so as to easily track progress

children <- 1:ncasem
for (icase in children) {
   hgt     <- hgtm[,icase]
   smoothList <-
      smooth.monotone(age, hgt, growfdPar, dbglev=0)
   Wfd     <- smoothList$Wfdobj
   beta    <- smoothList$beta
   Flist   <- smoothList$Flist
   iternum <- smoothList$iternum
   cvecm[,icase] <- Wfd$coefs
   betam[,icase] <- beta
   hgthat <- beta[1] + beta[2]*monfn(age, Wfd)
   RMSE   <- sqrt(mean((hgt - hgthat)^2*wgt)/mean(wgt))
   RMSEm[icase] <- RMSE
   cat(c(icase, iternum),paste("  ",round(Flist$f,4),
       "  ",round(RMSE, 4)))
}

# Females:

cvecf <- matrix(0, nbasis, ncasef)
betaf <- matrix(0, 2, ncasef)
RMSEf <- matrix(0, 1, ncasef)

children <- 1:ncasef
for (icase in children) {
   hgt    <- hgtf[,icase]
   smoothList <-
      smooth.monotone(age, hgt, growfdPar, dbglev=0)
   Wfd     <- smoothList$Wfd
   beta    <- smoothList$beta
   Flist   <- smoothList$Flist
   iternum <- smoothList$iternum
   cvecf[,icase] <- Wfd$coefs
   betaf[,icase] <- beta
   hgthat <- beta[1] + beta[2]*monfn(age, Wfd)
   RMSE   <- sqrt(mean((hgt - hgthat)^2*wgt)/mean(wgt))
   RMSEf[icase] <- RMSE
   cat(c(icase, iternum),paste("  ",round(Flist$f,4),
       "  ",round(RMSE, 4)))
}

#  -------------  plot the results  --------------------

#  blue:  monotone fit,  green:  non-monotone fit

#  Males:

par(mfrow=c(2,2),pty="s",ask=T)
children <- 1:ncasem
for (i in children) {
    #  curve values for monotone smooth
    Wfd  <- fd(cvecm[,i],wbasis)
    beta <- betam[,i]
    hgtmfit <- beta[1] + beta[2]*monfn(age, Wfd)
    hgtmhat <- beta[1] + beta[2]*monfn(agefine, Wfd)
    velmhat <- beta[2]*eval.monfd(agefine, Wfd, 1)
    accmhat <- beta[2]*eval.monfd(agefine, Wfd, 2)
    #  plot height data and fit 
    plot(age, hgtm[,i], ylim=c(60,200), 
         xlab="Years", ylab="", main=paste("Height for male",i))
    lines(agefine, hgtmhat, col=4)
    #  plot residuals
    resi <- hgtm[,i] - hgtmfit
    ind  <- resi >= -.7 & resi <= .7
    plot(age[ind], resi[ind], type="b", ylim=c(-.7,.7),
         xlab="Years", ylab="", main="Residuals")
    abline(h=0, lty=2)
    resiN <- hgtm[,i] - hgtmfitN[,i]
    indN  <- resiN >= -.7 & resiN <= .7
    points(age[indN], resiN[indN], col=3)
    lines(age[indN], resiN[indN], col=3)
    #  plot velocity
    ind <- velmhat >= 0 & velmhat <= 20
    plot(agefine[ind], velmhat[ind], type="l", ylim=c(0,20), col=4,
         xlab="Years", ylab="", main="Velocity")
    indN <- velmhatN[,i] >= 0 & velmhatN[,i] <= 20
    lines(agefine[indN], velmhatN[indN,i], col=3)
    #  plot acceleration
    ind <- accmhat >= -6 & accmhat <= 6
    plot(agefine[ind], accmhat[ind], type="l", ylim=c(-6,6), col=4,
         xlab="Years", ylab="", main="Acceleration")
    abline(h=0, lty=2)
    indN <- accmhatN[,i] >= -6 & accmhatN[,i] <= 6
    lines(agefine[indN], accmhatN[indN,i], col=3)
}

#  Females:

par(mfrow=c(2,2),pty="s",ask=T)
children <- 1:ncasef
for (i in children) {
    #  curve values for monotone smooth
    Wfd  <- fd(cvecf[,i],wbasis)
    beta <- betaf[,i]
    hgtffit <- beta[1] + beta[2]*monfn(age, Wfd)
    hgtfhat <- beta[1] + beta[2]*monfn(agefine, Wfd)
    velfhat <- beta[2]*eval.monfd(agefine, Wfd, 1)
    accfhat <- beta[2]*eval.monfd(agefine, Wfd, 2)
    #  plot height data and fit 
    plot(age, hgtf[,i], ylim=c(60,200), 
         xlab="Years", ylab="", main=paste("Height for female",i))
    lines(agefine, hgtfhat, col=4)
    #  plot residuals
    resi <- hgtf[,i] - hgtffit
    ind  <- resi >= -.7 & resi <= .7
    plot(age[ind], resi[ind], type="b", ylim=c(-.7,.7),
         xlab="Years", ylab="", main="Residuals")
    abline(h=0, lty=2)
    resiN <- hgtf[,i] - hgtffitN[,i]
    indN  <- resiN >= -.7 & resiN <= .7
    points(age[indN], resiN[indN], col=3)
    lines(age[indN], resiN[indN], col=3)
    #  plot velocity
    ind <- velfhat >= 0 & velfhat <= 20
    plot(agefine[ind], velfhat[ind], type="l", ylim=c(0,20), col=4,
         xlab="Years", ylab="", main="Velocity")
    indN <- velfhatN[,i] >= 0 & velfhatN[,i] <= 20
    lines(agefine[indN], velfhatN[indN,i], col=3)
    #  plot acceleration
    ind <- accfhat >= -6 & accfhat <= 6
    plot(agefine[ind], accfhat[ind], type="l", ylim=c(-6,6), col=4,
         xlab="Years", ylab="", main="Acceleration")
    abline(h=0, lty=2)
    indN <- accfhatN[,i] >= -6 & accfhatN[,i] <= 6
    lines(agefine[indN], accfhatN[indN,i], col=3)
}


