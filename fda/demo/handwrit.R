#  -----------------------------------------------------------------------
#                 Registered Handwriting Data
#  -----------------------------------------------------------------------
#
#                          Overview of the analyses
#
# These data are the X-Y coordinates of 20 replications of writing
# the script "fda".  The subject was Jim Ramsay.  Each replication
# is represented by 1401 coordinate values.  The scripts have been
# extensively pre-processed.  They have been adjusted to a common
# length that corresponds to 2.3 seconds or 2300 milliseconds, and
# they have already been registered so that important features in
# each script are aligned.
#
# This analysis is designed to illustrate techniques for working
# with functional data having rather high frequency variation and
# represented by thousands of data points per record.  Comments
# along the way explain the choices of analysis that were made.
#
# The final result of the analysis is a third order linear
# differential equation for each coordinate forced by a
# constant and by time.  The equations are able to reconstruct
# the scripts to a fairly high level of accuracy, and are also
# able to accommodate a substantial amount of the variation in
# the observed scripts across replications.  by contrast, a
# second order equation was found to be completely inadequate.
#
# An interesting suprise in the results is the role placed by
# a 120 millisecond cycle such that sharp features such as cusps
# correspond closely to this period.  This 110-120 msec cycle
# seems is usually seen in human movement data involving rapid
# movements, such as speech, juggling and so on.
#  -----------------------------------------------------------------------

#  Last modified 28 August 2012 by Jim Ramsay

#  Attach the FDA functions

library(fda)

#  Input the data.  These 20 records have already been
#  normalized to a common time interval of 2300 milliseconds
#  and have beeen also registered so that prominent features
#  occur at the same times across replications.
#  Time will be measured in milliseconds and space in meters.
#  The data will require a small amount of smoothing, since
#  an error of 0.5 mm is characteristic of the OPTOTRAK 3D
#  measurement system used to collect the data.

#  input the data

#temp <- array(scan("../data/fdareg.txt",0),c(20,2,1401))

#  set up a three-dimensional array

#fdaarray <- array(0, c(1401, 20, 2))
#fdaarray[,,1] <- t(temp[,1,])/1000
#fdaarray[,,2] <- t(temp[,2,])/1000
#imnames(fdaarray) <- list(NULL, NULL, c("X", "Y") )

#  Set up time values and range.
#  It is best to choose milliseconds as a time scale
#  in order to make the ratio of the time
#  unit to the inter-knot interval not too
#  far from one.  Otherwise, smoothing parameter values
#  may be extremely small or extremely large.

fdaarray = handwrit

fdatime  <- seq(0, 2300, len=1401)
fdarange <- c(0, 2300)

#  The basis functions will be B-splines, with a spline
#  placed at each knot.  One may question whether so many
#  basis functions are required, but this decision is found to
#  be  essential for stable derivative estimation up to the
#  third order at and near the boundaries.

#  Order 7 was used to get a smooth third derivative, which
#  requires penalizing the size of the 5th derivative, which
#  in turn requires an order of at least 7.
#  This implies norder + no. of interior knots = 1399 + 7 = 1406
#  basis functions.

nbasis   <- 1406
norder   <-    7
fdabasis <- create.bspline.basis(fdarange, nbasis, norder)

#  The smoothing parameter value 1e8 was chosen to obtain a
#  fitting error of about 0.5 mm, the known error level in
#  the OPTOTRACK equipment.

fdafd  <- fd(array(0, c(nbasis,20,2)), fdabasis)
lambda <- 1e8
fdaPar <- fdPar(fdafd, 5, lambda)

#  set up the functional data structure
#  (required 1.5 mins on a Lenovo X201 laptop)

smoothList <- smooth.basis(fdatime, fdaarray, fdaPar)
fdafd <- smoothList$fd
df    <- smoothList$df
gcv   <- smoothList$gcv
#  Add suitable names for the dimensions of the data.
fdafd$fdnames[[1]] <- "Milliseconds"
fdafd$fdnames[[2]] <- "Replications"
fdafd$fdnames[[3]] <- "Metres"
#  display degrees of freedom and total GCV criterion
df  #  about 115
totalgcv <- sum(gcv)
totalgcv
RMSgcv <- sqrt(totalgcv)*1000 # about 0.3 mm
RMSgcv

#  plot the fit to the data

par(mfrow=c(2,1),pty="m")
plotfit.fd(fdaarray, fdatime, fdafd)

#  plot all curves

par(mfrow=c(2,1),pty="m",ask=FALSE)
plot(fdafd)

#  compute values of curves and the values of the curve

fdameanfd  <- mean(fdafd)
fdamat     <- eval.fd(fdatime, fdafd)
fdameanmat <- apply(fdamat, c(1,3), mean)

#  Set up motor control clock cycle times at every
#  119 milliseconds.

cycle  <- seq(0,2300,119)
ncycle <- length(cycle)

#  evaluate curves at cycle times

fdamatcycle     <- eval.fd(cycle, fdafd)
fdameanmatcycle <- apply(fdamatcycle,c(1,3),mean)

#  Indices of cycle times corresponding to important features:
#  -- the cusp in "f",
#  -- the the cusp in "d",
#  -- the first cusp in "a",
#  -- the rest after the first cusp in "a", and
#  -- the second cusp in "a".
#  It is remarkable that these features correspond so closely
#  with clock cycle times!

featureindex   <- c(3, 5, 7, 10, 13, 16, 19)
fdafeature     <- fdamatcycle[featureindex,,]
fdameanfeature <- fdameanmatcycle[featureindex,]

#  Plot mean, including both sampling points and fit
#  Points at cycle times are plotted as blue circles, and
#  points at feature times are plotted as red circles.

par(mfrow=c(1,1), pty="s")
plot(fdameanmat[,1], fdameanmat[,2], type="l", lwd=2,
     xlab="Metres", ylab="Metres",
     xlim=c(-.040, .040), ylim=c(-.040, .040),
     main="Mean script")
points(fdameanmatcycle[-featureindex,1],
       fdameanmatcycle[-featureindex,2],        cex=1.2,col=2,lwd=4)
points(fdameanfeature[,1],  fdameanfeature[,2], cex=1.2,col=3,lwd=4)

#  Plot individual curves, including both sampling points and fit
#  also plot the mean curve in the background.
#  Note how closely individual curve features are tied to the
#  feature cycle times.

par(mfrow=c(1,1), pty="s",ask=TRUE)
for (i in 1:20) {
    plot(fdamat[,i,1], fdamat[,i,2], type="l", lwd=2,
         xlab="Metres", ylab="Metres",
         xlim=c(-.040, .040), ylim=c(-.040, .040),
         main=paste("Script",i))
    points(fdamatcycle[-featureindex,i,1],
          fdamatcycle[-featureindex,i,2],         cex=1.2,col=2,lwd=4)
    points(fdafeature[,i,1],    fdafeature[,i,2],  cex=1.2,col=3,lwd=4)
    lines( fdameanmat[,1],      fdameanmat[,2],      lty=4)
    points(fdameanmatcycle[-featureindex,1],
           fdameanmatcycle[-featureindex,2],         cex=1.2,col=2,lwd=4)
    points(fdameanfeature[,1],  fdameanfeature[,2],  cex=1.2,col=3,lwd=4)
}

#  Evaluate the three derivatives and their means

D1fdamat <- eval.fd(fdatime, fdafd, 1)
D2fdamat <- eval.fd(fdatime, fdafd, 2)
D3fdamat <- eval.fd(fdatime, fdafd, 3)

D1fdameanmat <- apply(D1fdamat, c(1,3), mean)
D2fdameanmat <- apply(D2fdamat, c(1,3), mean)
D3fdameanmat <- apply(D3fdamat, c(1,3), mean)

#  Plot the individual acceleration records.
#  In these plots, acceleration is displayed as
#  metres per second per second.

#  Cycle and feature times are plotted as vertical
#  dashed lines, un-featured cycle times as red
#  dotted lines, and cycle times of features as
#  heavier magenta solid lines.

par(mfrow=c(1,1), mar=c(5,5,4,2), pty="m",ask=TRUE)
for (i in 1:20) {
    matplot(fdatime, cbind(1e6*D2fdamat[,i,1],1e6*D2fdamat[,i,2]),
         type="l", lty=1, cex=1.2, col=c(2,4),
         xlim=c(0, 2300), ylim=c(-12, 12),
         xlab="Milliseconds", ylab="Meters/msec/msec",
         main=paste("Curve ",i))
    abline(h=0, lty=2)
    plotrange <- c(-12,12)
    for (k in 1:length(cycle)) abline(v=cycle[k], lty=2)
    for (j in 1:length(featureindex))
	     abline(v=cycle[featureindex[j]], lty=1)
    legend(1800,11.5, c("X", "Y"), lty=1, col=c(2,4))
}

#  Compute and plot the acceleration magnitudes,
#  also called the tangential accelerations.

D2mag  <- sqrt(D2fdamat[,,1]^2 + D2fdamat[,,2]^2)
D2magmean <- apply(D2mag,1,mean)

cexval <- 1.2
par(mfrow=c(1,1), mar=c(5,5,4,2)+cexval+2, pty="m",ask=FALSE)
matplot(fdatime, 1e6*D2mag, type="l", cex=1.2,
     xlab="Milliseconds", ylab="Metres/sec/sec",
     xlim=c(0,2300), ylim=c(0,12),
     main="Acceleration Magnitude")
plotrange <- c(0,12)
for (k in 1:length(cycle)) abline(v=cycle[k], lty=2)
for (j in 1:length(featureindex))
	abline(v=cycle[featureindex[j]], lty=1)

#  Plot the mean acceleration magnitude as well as
#  those for each curve.
#  Note the two rest cycles, one in "d" and one in "a"

par(mfrow=c(1,1), mar=c(5,5,4,2)+cexval+2, pty="m")
plot(fdatime, 1e6*D2magmean, type="l", cex=1.2,
     xlab="Milliseconds", ylab="Metres/sec/sec",
     xlim=c(0,2300), ylim=c(0,8),
     main="Mean acceleration Magnitude")
plotrange <- c(0,8)
for (k in 1:length(cycle)) abline(v=cycle[k], lty=2)
for (j in 1:length(featureindex))
	abline(v=cycle[featureindex[j]], lty=1)

#  Plot each individual acceleration magnitude, along
#  with the mean magnitude as a green dashed line

par(mfrow=c(1,1), mar=c(5,5,4,2), pty="m", ask=TRUE)
plotrange <- c(0,12)
for (i in 1:20) {
    plot(fdatime, 1e6*D2mag[,i], type="l", cex=1.2,
         xlim=c(0,2300), ylim=c(0,12),
         xlab="Milliseconds", ylab="Metres/sec/sec",
         main=paste("Script ",i))
    lines(fdatime,  1e6*D2magmean, lty=3)
    for (k in 1:length(cycle)) abline(v=cycle[k], lty=2)
    for (j in 1:length(featureindex))
	     abline(v=cycle[featureindex[j]], lty=1)
}

#  ------------------------------------------------------------
#                 Principal Differential Analysis
#  A third order equation forced by a constant and time is
#          estimated for X and Y coordinates separately.
#  Forcing with constant and time is required to allow for
#  an arbitrary origin and the left-right motion in X.
#  ------------------------------------------------------------

difeorder  <- 3  #  order of equation

#  set up the two forcing functions

ufdlist <- vector("list", 2)
#  constant forcing
constbasis <- create.constant.basis(fdarange)
constfd    <- fd(matrix(1,1,20), constbasis)
ufdlist[[1]] <- constfd
# time forcing
linbasis   <- create.monomial.basis(fdarange, 2)
lincoef    <- matrix(0,2,20)
lincoef[2,] <- 1
ufdlist[[2]] <- fd(lincoef, linbasis)

#  set up the corresponding weight functions

awtlist    <- vector("list", 2)
constfd    <- fd(1, constbasis)
constfdPar <- fdPar(constfd)
awtlist[[1]] <- constfdPar
awtlist[[2]] <- constfdPar

#  Define two basis systems for the derivative weight
#  functions.  One for a background analysis used as a
#  baseline and using constant weight functions, and
#  another using 125 basis functions that will be used
#  to generate the equaations.

#  Set the number of basis functions to 125,
#  found to maximize the RSQ measure, and corresponding
#  to DF above, the equivalent degrees of freedom in
#  the smooth.

wbasis125 <- create.bspline.basis(fdarange, 125)

#  ------------------------------------------------------------
#                   Analysis for coordinate X
#  ------------------------------------------------------------

#  Define the variable

fdafdX <- smooth.basis(fdatime, fdaarray[,,1], fdaPar)$fd

xfdlist <- vector("list", 1)
xfdlist[[1]] <- fdafdX

#  Set up the derivative weight functions

bwtlist <- vector("list", 3)
bfd     <- fd(matrix(0,1,1), constbasis)
bfdPar  <- fdPar(bfd, 1, 0)
bwtlist[[1]] <- bfdPar
bwtlist[[2]] <- bfdPar
bwtlist[[3]] <- bfdPar

#  carry out principal differential analysis
#  (this takes about 2 minutes on an IBM X31 notebook)

pdaList <- pda.fd(xfdlist, bwtlist, awtlist, ufdlist)

bestwtlist <- pdaList$bwtlist
aestwtlist <- pdaList$awtlist
resfdlist  <- pdaList$resfdlist

#  evaluate forcing functions

resfd  <- resfdlist[[1]]
resmat <- eval.fd(fdatime, resfd)

MSY <- mean(resmat^2)
MSY

#  Set the number of basis functions to 125,

bfd     <- fd(matrix(0,125,1), wbasis125)
bfdPar  <- fdPar(bfd, 1, 0)
bwtlist <- vector("list", 3)
bwtlist[[1]] <- bfdPar
bwtlist[[2]] <- bfdPar
bwtlist[[3]] <- bfdPar

#  carry out principal differential analysis
#  (this takes about 2 minutes on an Lenovo X201 laptop)

pdaList <- pda.fd(xfdlist, bwtlist, awtlist, ufdlist, difeorder)
bestwtlist <- pdaList$bwtlist
aestwtlist <- pdaList$awtlist
resfdlist  <- pdaList$resfdlist

#  evaluate forcing functions

resfd  <- resfdlist[[1]]
resmat <- eval.fd(fdatime, resfd)

#  compute a squared multiple correlation measure of fit
#  MSY = mean(mean(resmat^2)) #  Used only with constant basis
#  Uncomment this line when the constant basis is used, and
#  comment it out otherwise.

MSE <- mean(resmat^2)
RSQ <- (MSY-MSE)/MSY
RSQ

#  Plot the weight functions

par(mfrow=c(3,1), ask=FALSE)
for (j in 1:3) {
	 betafdPar <- bestwtlist[[j]]
    plot(betafdPar$fd, cex=1, ylab=paste("Weight function ",j-1))
}

#  Plot the second derivative weight, defining the period
#  of a harmonic oscillator.

b2fdParX <- bestwtlist[[2]]
b2fdX    <- b2fdParX$fd
b2vecX   <- eval.fd(fdatime, b2fdX)
b2meanX  <- mean(b2vecX)

par(mfrow=c(1,1), pty="m")
plot(fdatime, b2vecX, type="l", cex=1.2,
     xlim=c(0, 2300), ylim=c(0, 6e-3))
abline(h=b2meanX, lty=3)
plotrange <- c(0,6e-3)
for (k in 1:length(cycle)) abline(v=cycle[k], lty=2)
for (j in 1:length(featureindex))
	abline(v=cycle[featureindex[j]], lty=1)

#  display coefficients for forcing weight functions

aestwtlist[[1]]$fd$coefs
aestwtlist[[2]]$fd$coefs

#  plot all forcing functions

par(mfrow=c(1,1), pty="m")
matplot(fdatime, 1e9*resmat, type="l", cex=1.2,
     xlim=c(0,2300), ylim=c(-200,200),
     xlab="Milliseconds", ylab="Meters/sec/sec/sec")
lines(fdatime, 1e9*D3fdameanmat[,1], lty=2)

#  plot the mean forcing function along with third deriv.

resmeanfd  <- mean(resfd)
resmeanvec <- eval.fd(fdatime, resmeanfd)

par(mfrow=c(1,1), pty="m")
plot(fdatime, 1e9*resmeanvec, type="l", cex=1.2, col=2,
     xlim=c(0,2300), ylim=c(-200,200),
     xlab="Milliseconds", ylab="Meters/sec/sec/sec")
lines(fdatime, 1e9*D3fdameanmat[,1], lty=2, col=3)

# Define a functional data object for the
#  three derivative weight functions

wcoef1 <- bestwtlist[[1]]$fd$coefs
wcoef2 <- bestwtlist[[2]]$fd$coefs
wcoef3 <- bestwtlist[[3]]$fd$coefs

wcoef  <- cbind(wcoef1, wcoef2, wcoef3)
wfd    <- fd(wcoef,wbasis125)

#  Set up a linear differential operator.
#  This isn"t used in these analyses.

fdaLfd <- Lfd(difeorder, fd2list(wfd))

ystart <- matrix(0,3,3)
ystart[1,1] <-   fdameanmat[1,1]
ystart[2,2] <- D1fdameanmat[1,1]
ystart[3,3] <- D2fdameanmat[1,1]

#  solve the equation

#  (Solving the equations with the following error tolerance takes 20 minutes
#   on an IBM X31 notebook, and about 5 seconds in Matlab.  Swtich to Matlab
#   if you have a lot of this kind of work to do.)

EPSval = 1e-4
odeList <- odesolv(bestwtlist, ystart, EPS=EPSval, MAXSTP=1e6)
tX <- odeList[[1]]
yX <- odeList[[2]]

#  plot the three solutions

par(mfrow=c(3,1), pty="m")

pltrng <- c(min(yX[1,,]), max(yX[1,,]))
matplot(tX, t(yX[1,,]), type="l", lty=1, ylim=pltrng, main="Function")
abline(h=0, lty=2)

pltrng <- c(min(yX[2,,]), max(yX[2,,]))
matplot(tX, t(yX[2,,]), type="l", lty=1, ylim=pltrng, main="Derivative")
abline(h=0, lty=2)

pltrng <- c(min(yX[3,,]), max(yX[3,,]))
matplot(tX, t(yX[3,,]), type="l", lty=1, ylim=pltrng, main="Derivative")
abline(h=0, lty=2)

#  set up curve and derivative values

umatx <- matrix(0,length(fdatime),3)
umatx[,1] <- approx(tX, t(yX[1,1,]), fdatime)$y
umatx[,2] <- approx(tX, t(yX[1,2,]), fdatime)$y
umatx[,3] <- approx(tX, t(yX[1,2,]), fdatime)$y

Dumatx <- matrix(0,length(fdatime),3)
Dumatx[,1] <- approx(tX, t(yX[2,1,]), fdatime)$y
Dumatx[,2] <- approx(tX, t(yX[2,2,]), fdatime)$y
Dumatx[,3] <- approx(tX, t(yX[2,3,]), fdatime)$y

D2umatx <- matrix(0,length(fdatime),3)
D2umatx[,1] <- approx(tX, t(yX[3,1,]), fdatime)$y
D2umatx[,2] <- approx(tX, t(yX[3,2,]), fdatime)$y
D2umatx[,3] <- approx(tX, t(yX[3,3,]), fdatime)$y

#  plot fit to each curve

par(mfrow=c(1,1), ask=TRUE, pty="m")
index  <- 1:20
fdamat <- eval.fd(fdatime, fdafd)
zmat   <- cbind(fdatime-1150,umatx)
for (i in index) {
    xhat <- fdamat[,i,1] - lsfit(zmat, fdamat[,i,1])$residual
    matplot(fdatime, cbind(xhat, fdamat[,i,1]),
            type="l", lty=c(1,3), cex=1.2,
            xlim=c(0, 2300), ylim=c(-0.04, 0.04),
            main=paste("X curve ",i))
}

#  ------------------------------------------------------------
#                   Analysis for coordinate Y
#  ------------------------------------------------------------

#  Define the variable

fdafdY <- smooth.basis(fdatime, fdaarray[,,2], fdaPar)$fd

yfdlist <- vector("list", 1)
yfdlist[[1]] <- fdafdY

#  Set up the derivative weight functions

bwtlist <- vector("list", 3)
bfd     <- fd(matrix(0,1,1), constbasis)
bfdPar  <- fdPar(bfd, 1, 0)
bwtlist[[1]] <- bfdPar
bwtlist[[2]] <- bfdPar
bwtlist[[3]] <- bfdPar

#  carry out principal differential analysis
#  (this takes about 2 minutes on an IBM X31 notebook)

pdaList <- pda.fd(yfdlist, bwtlist, awtlist, ufdlist, difeorder)
bestwtlist <- pdaList$bwtlist
aestwtlist <- pdaList$awtlist
resfdlist  <- pdaList$resfdlist

#  evaluate forcing functions

resfd  <- resfdlist[[1]]
resmat <- eval.fd(fdatime, resfd)

MSY <- mean(resmat^2)
MSY

#  Set the number of basis functions to 125,

bfd     <- fd(matrix(0,125,1), wbasis125)
bfdPar  <- fdPar(bfd, 1, 0)
bwtlist <- vector("list", 3)
bwtlist[[1]] <- bfdPar
bwtlist[[2]] <- bfdPar
bwtlist[[3]] <- bfdPar

#  carry out principal differential analysis
#  (this takes about 2 minutes on an IBM X31 Notebook)

pdaList <- pda.fd(yfdlist, bwtlist, awtlist, ufdlist, difeorder)
bestwtlist <- pdaList$bwtlist
aestwtlist <- pdaList$awtlist
resfdlist  <- pdaList$resfdlist

#  evaluate forcing functions

resfd  <- resfdlist[[1]]
resmat <- eval.fd(fdatime, resfd)

#  compute a squared multiple correlation measure of fit
#  MSY = mean(mean(resmat^2)) #  Used only with constant basis
#  Uncomment this line when the constant basis is used, and
#  comment it out otherwise.

MSE <- mean(resmat^2)
RSQ <- (MSY-MSE)/MSY
RSQ

#  Plot the weight functions

par(mfrow=c(3,1),ask=FALSE)
for (j in 1:3) {
	 betafdPar <- bestwtlist[[j]]
    plot(betafdPar$fd, cex=1, ylab=paste("Weight function ",j-1))
}

#  Plot the second derivative weight, defining the period
#  of a harmonic oscillator.

b2fdParY <- bestwtlist[[2]]
b2fdY    <- b2fdParY$fd
b2vecY   <- eval.fd(fdatime, b2fdY)
b2meanY  <- mean(b2vecY)

par(mfrow=c(1,1), pty="m", ask=FALSE)
plot(fdatime, b2vecY, type="l", cex=1.2,
     xlim=c(0, 2300), ylim=c(0, 6e-3))
abline(h=b2meanY, lty=3)
plotrange <- c(0,6e-3)
for (k in 1:length(cycle)) abline(v=cycle[k], lty=2)
for (j in 1:length(featureindex))
	abline(v=cycle[featureindex[j]], lty=1)

#  display coefficients for forcing weight functions

aestwtlist[[1]]$fd$coefs
aestwtlist[[2]]$fd$coefs

#  plot all forcing functions

par(mfrow=c(1,1), ask=FALSE, pty="m")
matplot(fdatime, 1e9*resmat, type="l", cex=1.2,
     xlim=c(0,2300), ylim=c(-200,200),
     xlab="Milliseconds", ylab="Meters/sec/sec/sec")
lines(fdatime, 1e9*D3fdameanmat[,1], lty=2)

#  plot the mean forcing function along with third deriv.

resmeanfd  <- mean(resfd)
resmeanvec <- eval.fd(fdatime, resmeanfd)

par(mfrow=c(1,1), ask=FALSE, pty="m")
plot(fdatime, 1e9*resmeanvec, type="l", cex=1.2, col=2,
     xlim=c(0,2300), ylim=c(-200,200),
     xlab="Milliseconds", ylab="Meters/sec/sec/sec")
lines(fdatime, 1e9*D3fdameanmat[,1], lty=2, col=3)

# Define a functional data object for the
#  three derivative weight functions

wcoef1 <- bestwtlist[[1]]$fd$coefs
wcoef2 <- bestwtlist[[2]]$fd$coefs
wcoef3 <- bestwtlist[[3]]$fd$coefs

wcoef  <- cbind(wcoef1, wcoef2, wcoef3)
wfd    <- fd(wcoef,wbasis125)

#  Set up a linear differential operator.
#  This isn"t used in these analyses.

fdaLfd <- Lfd(difeorder, fd2list(wfd))

ystart <- matrix(0,3,3)
ystart[1,1] <-   fdameanmat[1,2]
ystart[2,2] <- D1fdameanmat[1,2]
ystart[3,3] <- D2fdameanmat[1,2]

#  solve the equation

#  (Solving the equations with the following error tolerance takes 20 minutes
#   on an IBM X31 notebook, and about 5 seconds in Matlab.  Swtich to Matlab
#   if you have a lot of this kind of work to do.)

EPSval = 1e-4
odeList <- odesolv(bestwtlist, ystart, EPS=EPSval, MAXSTP=1e6)
tY <- odeList[[1]]
yY <- odeList[[2]]

#  plot the three solutions

par(mfrow=c(3,1), ask=FALSE, pty="m")

pltrng <- c(min(yY[1,,]), max(yY[1,,]))
matplot(tY, t(yY[1,,]), type="l", lty=1, ylim=pltrng, main="Function")
abline(h=0, lty=2)

pltrng <- c(min(yY[2,,]), max(yY[2,,]))
matplot(tY, t(yY[2,,]), type="l", lty=1, ylim=pltrng, main="Derivative")
abline(h=0, lty=2)

pltrng <- c(min(yY[3,,]), max(yY[3,,]))
matplot(tY, t(yY[3,,]), type="l", lty=1, ylim=pltrng, main="Derivative")
abline(h=0, lty=2)

#  set up curve and derivative values

umaty <- matrix(0,length(fdatime),3)
umaty[,1] <- approx(tY, t(yY[1,1,]), fdatime)$y
umaty[,2] <- approx(tY, t(yY[1,2,]), fdatime)$y
umaty[,3] <- approx(tY, t(yY[1,2,]), fdatime)$y

Dumaty <- matrix(0,length(fdatime),3)
Dumaty[,1] <- approx(tY, t(yY[2,1,]), fdatime)$y
Dumaty[,2] <- approx(tY, t(yY[2,2,]), fdatime)$y
Dumaty[,3] <- approx(tY, t(yY[2,3,]), fdatime)$y

D2umaty <- matrix(0,length(fdatime),3)
D2umaty[,1] <- approx(tY, t(yY[3,1,]), fdatime)$y
D2umaty[,2] <- approx(tY, t(yY[3,2,]), fdatime)$y
D2umaty[,3] <- approx(tY, t(yY[3,3,]), fdatime)$y

#  plot fit to each curve

par(mfrow=c(1,1), ask=TRUE, pty="m")
index  <- 1:20
fdamat <- array(0,c(1401,20,2))
fdamat[,,1] <- eval.fd(fdatime, fdafdX)
fdamat[,,2] <- eval.fd(fdatime, fdafdY)
zmat   <- cbind(fdatime-1150,umaty)
for (i in index) {
    yhat <- fdamat[,i,2] - lsfit(zmat, fdamat[,i,2])$residual
    matplot(fdatime, cbind(yhat, fdamat[,i,2]),
            type="l", lty=c(1,3), cex=1.2,
            xlim=c(0, 2300), ylim=c(-0.04, 0.04),
            main=paste("Y curve ",i))
}

#  plot the two weight functions for the second derivative

par(mfrow=c(1,1), mar=c(5,5,4,2)+cexval+2, pty="m")
matplot(fdatime, cbind(b2vecX, b2vecY), type="l", lty=1, col=c(2,4),
        xlim=c(0, 2300), ylim=c(0, 6e-3))
abline(h=b2meanX, lty=3, col=2)
abline(h=b2meanY, lty=3, col=4)
plotrange <- c(0,6e-3)
for (k in 1:length(cycle)) abline(v=cycle[k], lty=2)
for (j in 1:length(featureindex))
	abline(v=cycle[featureindex[j]], lty=1)


