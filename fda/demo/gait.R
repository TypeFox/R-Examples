#  --------------------------------------------------------------------
#                        Gait data
#  --------------------------------------------------------------------

#  --------------------------------------------------------------------
#
#                          Overview of the analyses
#
#  The gait data were chosen for these sample analyses because they are
#  bivariate:  consisting of both hip and knee angles observed over a
#  gait cycle for 39 children.  The bivariate nature of the data implies
#  certain displays and analyses that are not usually considered, and
#  especially the use of canonical correlation analysis.
#
#  As with the daily weather data, the harmonic acceleration roughness
#  penalty is used throughout since the data are periodic with a strong
#  sinusoidal component of variation.
#
#  After setting up the data, smoothing the data using GCV (generalized
#  cross-validation) to select a smoothing parameter, and displaying
#  various descriptive results, the data are subjected to a principal
#  components analysis, followed by a canonical correlation analysis of
#  thejoint variation of hip and knee angle, and finally a registration
#  of the curves.  The registration is included here especially because
#  the registering of periodic data requires the estimation of a phase
#  shift constant for each curve in addition to possible nonlinear
#  transformations of time.
#
#  --------------------------------------------------------------------

#  Last modified 10 November 2010 by Jim Ramsay

#  attach the FDA functions

library(fda)

#  Set up the argument values: equally spaced over circle of
#  circumference 20.  Earlier  analyses of the gait data used time
#  values over [0,1], but led to singularity problems in the use of
#  function fRegress.  In general, it is better use a time interval
#  that assigns roughly one time unit to each inter-knot interval.

(gaittime <- as.numeric(dimnames(gait)[[1]])*20)
gaitrange <- c(0,20)

#  display ranges of gait for each variable

apply(gait, 3, range)

# -----------  set up the harmonic acceleration operator  ----------

harmaccelLfd <- vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval=gaitrange)

#  Set up basis for representing gait data.  The basis is saturated
#  since there are 20 data points per curve, and this set up defines
#  21 basis functions.  Recall that a fourier basis has an odd number
#  of basis functions.

gaitbasis <- create.fourier.basis(gaitrange, nbasis=21)

#  -------------------------------------------------------------------
#                 Choose level of smoothing using
#          the generalized cross-validation criterion
#  -------------------------------------------------------------------

#  set up range of smoothing parameters in log_10 units

gaitLoglam <- seq(-4,0,0.25)
nglam   <- length(gaitLoglam)

gaitSmoothStats <- array(NA, dim=c(nglam, 3),
      dimnames=list(gaitLoglam, c("log10.lambda", "df", "gcv") ) )
gaitSmoothStats[, 1] <- gaitLoglam

#  loop through smoothing parameters

for (ilam in 1:nglam) {
  gaitSmooth <- smooth.basisPar(gaittime, gait, gaitbasis,
                   Lfdobj=harmaccelLfd, lambda=10^gaitLoglam[ilam])
  gaitSmoothStats[ilam, "df"]  <- gaitSmooth$df
  gaitSmoothStats[ilam, "gcv"] <- sum(gaitSmooth$gcv)
  # note: gcv is a matrix in this case
}

#  display and plot GCV criterion and degrees of freedom

gaitSmoothStats
plot(gaitSmoothStats[, 1], gaitSmoothStats[, 3])

#  set up plotting arrangements for one and two panel displays
#  allowing for larger fonts

op <- par(mfrow=c(2,1))
plot(gaitLoglam, gaitSmoothStats[, "gcv"], type="b",
     xlab="Log_10 lambda", ylab="GCV Criterion",
     main="Gait Smoothing", log="y")

plot(gaitLoglam, gaitSmoothStats[, "df"], type="b",
     xlab="Log_10 lambda", ylab="Degrees of freedom",
     main="Gait Smoothing")
par(op)

# With gaittime <- (1:20)/21,
#    GCV is minimized with lambda = 10^(-2).

gaitfd <- smooth.basisPar(gaittime, gait,
       gaitbasis, Lfdobj=harmaccelLfd, lambda=1e-2)$fd
names(gaitfd$fdnames) <- c("Normalized time", "Child", "Angle")
gaitfd$fdnames[[3]] <- c("Hip", "Knee")

str(gaitfd)

#  --------  plot curves and their first derivatives  ----------------

#par(mfrow=c(1,2), mar=c(3,4,2,1), pty="s")
op <- par(mfrow=c(2,1))
plot(gaitfd, cex=1.2)
par(op)

#  plot each pair of curves interactively

plotfit.fd(gait, gaittime, gaitfd, cex=1.2, ask=FALSE)

#  plot the residuals, sorting cases by residual sum of squares
#  this produces 39 plots for each of knee and hip angle

plotfit.fd(gait, gaittime, gaitfd, residual=TRUE, sort=TRUE, cex=1.2)

#  plot first derivative of all curves

op <- par(mfrow=c(2,1))
plot(gaitfd, Lfdobj=1)
par(op)

#  -----------------------------------------------------------------
#            Display the mean, variance and covariance functions
#  -----------------------------------------------------------------

#  ------------  compute the mean functions  --------------------

gaitmeanfd <- mean.fd(gaitfd)

#  plot these functions and their first two derivatives

op <- par(mfcol=2:3)
plot(gaitmeanfd)
plot(gaitmeanfd, Lfdobj=1)
plot(gaitmeanfd, Lfdobj=2)
par(op)

#  --------------  Compute the variance functions  -------------

gaitvarbifd <- var.fd(gaitfd)
str(gaitvarbifd)

gaitvararray <- eval.bifd(gaittime, gaittime, gaitvarbifd)

#  plot variance and covariance functions as contours

filled.contour(gaittime, gaittime, gaitvararray[,,1,1], cex=1.2)
title("Knee - Knee")

filled.contour(gaittime, gaittime, gaitvararray[,,1,2], cex=1.2)
title("Knee - Hip")

filled.contour(gaittime, gaittime, gaitvararray[,,1,3], cex=1.2)
title("Hip - Hip")

#  plot variance and covariance functions as surfaces

persp(gaittime, gaittime, gaitvararray[,,1,1], cex=1.2)
title("Knee - Knee")

persp(gaittime, gaittime, gaitvararray[,,1,2], cex=1.2)
title("Knee - Hip")

persp(gaittime, gaittime, gaitvararray[,,1,3], cex=1.2)
title("Hip - Hip")

#  plot correlation functions as contours

gaitCorArray <- cor.fd(gaittime, gaitfd)

quantile(gaitCorArray)

contour(gaittime, gaittime, gaitCorArray[,,1,1], cex=1.2)
title("Knee - Knee")

contour(gaittime, gaittime, gaitCorArray[,,1,2], cex=1.2)
title("Knee - Hip")

contour(gaittime, gaittime, gaitCorArray[,,1,3], cex=1.2)
title("Hip - Hip")

#  --------------------------------------------------------------
#            Principal components analysis
#  --------------------------------------------------------------

#  do the PCA with varimax rotation

# Smooth with lambda as determined above

gaitfdPar  <- fdPar(gaitbasis, harmaccelLfd, lambda=1e-2)

gaitpca.fd <- pca.fd(gaitfd, nharm=4, gaitfdPar)

gaitpca.fd <- varmx.pca.fd(gaitpca.fd)

#  plot harmonics using cycle plots

op <- par(mfrow=c(2,2))
plot.pca.fd(gaitpca.fd, cycle=TRUE)
par(op)

#  compute proportions of variance associated with each angle

#  compute the harmonic scores associated with each angle

gaitscores = gaitpca.fd$scores

#  compute the values of the harmonics at time values for each angle

gaitharmmat = eval.fd(gaittime, gaitpca.fd$harmonics)
hipharmmat  = gaitharmmat[,,1]
kneeharmmat = gaitharmmat[,,2]

#  we need the values of the two mean functions also

gaitmeanvec = eval.fd(gaittime, gaitmeanfd)
hipmeanvec  = gaitmeanvec[,,1]
kneemeanvec = gaitmeanvec[,,2]

#  the values of the smooths of each angle less each mean function

gaitsmtharray = eval.fd(gaittime, gaitfd)
hipresmat  = gaitsmtharray[,,1] - outer( hipmeanvec,rep(1,39))
kneeresmat = gaitsmtharray[,,2] - outer(kneemeanvec,rep(1,39))

#  the variances of the residuals of the smooth angles from their means

hipvar  = mean( hipresmat^2)
kneevar = mean(kneeresmat^2)

print(paste("Variances of fits by the means:",
            round(c(hipvar, kneevar),1)))

#  compute the fits to the residual from the mean achieved by the PCA

hipfitarray  = array(NA, c(nrow(hipharmmat ),nrow(gaitscores),ncol(gaitscores)))
kneefitarray = array(NA, c(nrow(kneeharmmat),nrow(gaitscores),ncol(gaitscores)))
for (isc in 1:2) {
               hipfitarray[,,isc]  = hipharmmat  %*% t(gaitscores[,,isc])
               kneefitarray[,,isc] = kneeharmmat %*% t(gaitscores[,,isc])
}

#  compute the variances of the PCA fits

hipfitvar  = c()
kneefitvar = c()
for (isc in 1:2) {
            hipfitvar  = c(hipfitvar,  mean( hipfitarray[,,isc]^2))
            kneefitvar = c(kneefitvar, mean(kneefitarray[,,isc]^2))
}

#  compute percentages relative to the total PCA fit variance
#  these percentages will add to 100

hippropvar1 = c()
kneepropvar1 = c()
for (isc in 2) {
            hippropvar1  = c(hippropvar1,  hipfitvar[isc]/(hipfitvar[isc] +
                                           kneefitvar[isc]))
            kneepropvar1 = c(kneepropvar1, kneefitvar[isc]/(hipfitvar[isc]+
                                           kneefitvar[isc]))
}

print(paste("Percentages of fits for the PCA:",
            round(100*c(hippropvar1, kneepropvar1),1)))

#  compute percentages relative to the total mean fit variance
#  these percentages will add to the total percentage of fit
#  accounted for by the pca, which will typically be less than 100

hippropvar2  = c()
kneepropvar2 = c()
for (isc in 1:2) {
            hippropvar2  = c(hippropvar2,  hipfitvar[isc] /(hipvar+kneevar))
            kneepropvar2 = c(kneepropvar2, kneefitvar[isc]/(hipvar+kneevar))
}

print((paste("Percentages of fits for the PCA:",
             round(100*c(hippropvar2, kneepropvar2),1))))

#  --------------------------------------------------------------
#           Canonical correlation analysis
#  --------------------------------------------------------------

hipfd  <- gaitfd[,1]
kneefd <- gaitfd[,2]

hipfdPar  <- fdPar(hipfd,  harmaccelLfd, 1e2)
kneefdPar <- fdPar(kneefd, harmaccelLfd, 1e2)

ccafd    <- cca.fd(hipfd, kneefd, ncan=3, hipfdPar, kneefdPar)

#  plot the canonical weight functions

op <- par(mfrow=c(2,1), mar=c(3,4,2,1), pty="m")
plot.cca.fd(ccafd, cex=1.2)
par(op)

#  display the canonical correlations

round(ccafd$ccacorr[1:6],3)
plot(1:6, ccafd$ccacorr[1:6], type="b")


#  --------------------------------------------------------------
#         Register the angular acceleration of the gait data
#  --------------------------------------------------------------

#  compute the acceleration and mean acceleration

D2gaitfd      <- deriv.fd(gaitfd,2)
names(D2gaitfd$fdnames)[[3]] <- "Angular acceleration"
D2gaitfd$fdnames[[3]] <- c("Hip", "Knee")
D2gaitmeanfd  <- mean.fd(D2gaitfd)
names(D2gaitmeanfd$fdnames)[[3]] <- "Mean angular acceleration"
D2gaitmeanfd$fdnames[[3]] <- c("Hip", "Knee")

#  set up basis for warping function

nwbasis   <- 7
wbasis    <- create.bspline.basis(gaitrange,nwbasis,3)
Warpfd    <- fd(matrix(0,nwbasis,5),wbasis)
WarpfdPar <- fdPar(Warpfd)

#  register the functions

gaitreglist <- register.fd(D2gaitmeanfd, D2gaitfd[1:5,], WarpfdPar, periodic=TRUE)

plotreg.fd(gaitreglist)

#  display horizonal shift values
print(round(gaitreglist$shift,1))

#  histogram of horizontal shift values
par(mfrow=c(1,1))
hist(gaitreglist$shift,xlab="Normalized time")

#  --------------------------------------------------------------
#              Predict knee angle from hip angle
#             for angle and angular acceleration
#  --------------------------------------------------------------

#  set up the data

hipfd  <- gaitfd[,1]
kneefd <- gaitfd[,2]
ncurve <- dim(kneefd$coefs)[2]

kneemeanfd <- mean(kneefd)

#  define the functional parameter object for regression functions

betafdPar <- fdPar(gaitbasis, harmaccelLfd)
betalist  <- list(betafdPar,betafdPar)

#  ----------  predict knee angle from hip angle --------

conbasis <- create.constant.basis(c(0,20))
constfd  <- fd(matrix(1,1,ncurve), conbasis)

#  set up the list of covariate objects

xfdlist  <- list(constfd, hipfd)

#  fit the current functional linear model

fRegressout <- fRegress(kneefd, xfdlist, betalist)

#  set up and plot the fit functions and the regression functions

kneehatfd   <- fRegressout$yhatfd
betaestlist <- fRegressout$betaestlist

alphafd   <- betaestlist[[1]]$fd
hipbetafd <- betaestlist[[2]]$fd

op <- par(mfrow=c(2,1), ask=FALSE)
plot(alphafd,   ylab="Intercept")
plot(hipbetafd, ylab="Hip coefficient")
par(op)

#  compute and plot squared multiple correlation function

gaitfine    <- seq(0,20,len=101)
kneemat     <- eval.fd(gaitfine, kneefd)
kneehatmat  <- predict(kneehatfd, gaitfine)
kneemeanvec <- as.vector(eval.fd(gaitfine, kneemeanfd))

SSE0 <- apply((kneemat - outer(kneemeanvec, rep(1,ncurve)))^2, 1, sum)
SSE1 <- apply((kneemat - kneehatmat)^2, 1, sum)
Rsqr <- (SSE0-SSE1)/SSE0

op <- par(mfrow=c(1,1),ask=FALSE)
plot(gaitfine, Rsqr, type="l", ylim=c(0,0.4))

#  for each case plot the function being fit, the fit,
#                     and the mean function

op <- par(mfrow=c(1,1),ask=TRUE)
for (i in 1:ncurve) {
  plot( gaitfine, kneemat[,i], type="l", lty=1, col=4, ylim=c(0,80))
  lines(gaitfine, kneemeanvec,           lty=2, col=2)
  lines(gaitfine, kneehatmat[,i],        lty=3, col=4)
  title(paste("Case",i))
}
par(op)

#  ----------  predict knee acceleration from hip acceleration --------

D2kneefd     <- deriv(kneefd, 2)
D2hipfd      <- deriv(hipfd, 2)
D2kneemeanfd <- mean(D2kneefd)

#  set up the list of covariate objects

D2xfdlist  <- list(constfd,D2hipfd)

#  fit the current functional linear model

D2fRegressout <- fRegress(D2kneefd, D2xfdlist, betalist)

#  set up and plot the fit functions and the regression functions

D2kneehatfd   <- D2fRegressout$yhatfd
D2betaestlist <- D2fRegressout$betaestlist

D2alphafd   <- D2betaestlist[[1]]$fd
D2hipbetafd <- D2betaestlist[[2]]$fd

op <- par(mfrow=c(2,1), ask=FALSE)
plot(D2alphafd,   ylab="D2Intercept")
plot(D2hipbetafd, ylab="D2Hip coefficient")
par(op)

#  compute and plot squared multiple correlation function

D2kneemat     <- eval.fd(gaitfine, D2kneefd)
D2kneehatmat  <- predict(D2kneehatfd, gaitfine)
D2kneemeanvec <- as.vector(eval.fd(gaitfine, D2kneemeanfd))

D2SSE0 <- apply((D2kneemat - outer(D2kneemeanvec, rep(1,ncurve)))^2, 1, sum)
D2SSE1 <- apply((D2kneemat - D2kneehatmat)^2, 1, sum)
D2Rsqr <- (D2SSE0-D2SSE1)/D2SSE0

par(mfrow=c(1,1),ask=FALSE)
plot(gaitfine, D2Rsqr, type="l", ylim=c(0,0.5))

#  for each case plot the function being fit, the fit, and the mean function

op <- par(mfrow=c(1,1),ask=TRUE)
for (i in 1:ncurve) {
  plot( gaitfine, D2kneemat[,i], type="l", lty=1, col=4, ylim=c(-20,20))
  lines(gaitfine, D2kneemeanvec,           lty=2, col=2)
  lines(gaitfine, D2kneehatmat[,i],        lty=3, col=4)
  lines(c(0,20), c(0,0), lty=2, col=2)
  title(paste("Case",i))
}
par(op)





