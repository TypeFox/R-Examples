
#  -----------------------------------------------------------------------
#                     Canadian Daily Weather Data
#  -----------------------------------------------------------------------

#  -----------------------------------------------------------------------
#
#                          Overview of the analyses
#
#  These analyses of the daily temperature and precipitation data for 35
#  Canadian weather stations are designed to be a tour of much of the
#  R and S-PLUS software.  The emphasis is on illustrating useful
#  techniques in a way that will help set up your own application.
#  Sophistication, computational efficiency and compression of code have
#  been ignored as considerations in favor of making the analyses as
#  transparent as possible.
#  The sections in this file are as follows:

#  As always, the first step is to inspect the data and to smooth the
#  data with the appropriate level of smoothing for the analyses being
#  considered.

#  The idea of a customized or "smart" roughness penalty
#  is introduced right away in the form of harmonic acceleration, which
#  is especially important for periodic data such as these with a
#  variation that is dominated by a sinusoid plus a constant signal, or
#  shifted harmonic variation.  However, in order to keep the analyses
#  simple and to economize on computational effort, we often compromise
#  the principle supported in the FDA book of using a saturated basis
#  capable of interpolating the data combined with a roughness penalty,
#  and instead opt for a Fourier series basis system with 65 basis
#  functions and no roughness penalty.

#  Nevertheless, there is a smoothing section below that uses 365
#  Fourier basis functions combined with a harmonic acceleration penalty
#  where we estimate the smoothing parameter by minimizing the
#  GCV or generalized cross-validation parameter.

#  Smoothing is followed by the display of various descriptive statistics,
#  including mean and standard deviation functions,  and covariance and
#  correlation surfaces.

#  A short section illustrates the principal components analysis of
#  temperature, and you may want to repeate these analyses for
#  precipitation.

#  The functional linear model is used in various ways in the following
#  sections:
#  1.  Functional analysis of variance is used to show climate zone
#      effects for temperature and precipitation.
#  2.  The concurrent functional linear model is illustrated by
#      predicting log precipitation from climate zone and a functional
#      covariate constructed by removing climate effects from temperature.
#  3.  Log annual precipitation, a scalar dependent variable, is fit
#      by using temperature as a functional covariate.  Harmonic
#      acceleration roughness in the regression coefficient function
#      is penalized.
#  4.  The full log precipitation function is fit by the regressing
#      on the full temperature profile, and various levels of smoothing
#      are used to show the the effects of smoothing over both arguments
#      of the regression coefficient function.

#  The final section illustrates the smoothing of precipitation by a
#  curve that is constrained to be strictly positive, as is required
#  for a variable like precipitation.
#  -----------------------------------------------------------------------

#  Last modified 5 November 2008 by Jim Ramsay
#  Previously modified 2 March 2007 by Spencer Graves

#  ------------------------  input the data  -----------------------

#  set up the times of observation at noon

#  -------------  set up fourier basis  ---------------------------
#  Here it was decided that 65 basis functions captured enough of
#  the detail in the temperature data: about one basis function
#  per week.

#  The use of only 65 basis functions instead of 365
#  automatically generates some smoothing.

#  However, see below for smoothing with a saturated
#  basis (365 basis functions) where smoothing is defined by the
#  GCV criterion.

daybasis65 <- create.fourier.basis(rangeval=c(0, 365), nbasis=65)

#  -----------  set up the harmonic acceleration operator  ----------

harmaccelLfd365 <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))

#  ---------  create fd objects for temp. and prec. ---------------

# First check the distribution
qqnorm(CanadianWeather$dailyAv[,,"Temperature.C"], datax=TRUE)
# Consistent with a strong annual cycle
# plus weaker normal noise

daytempfd <- with(CanadianWeather, smooth.basis(day.5,
       dailyAv[,,"Temperature.C"],
       daybasis65, fdnames=list("Day", "Station", "Deg C"))$fd )
plot(daytempfd, axes=FALSE)
axisIntervals(1)
axis(2)

# Check precipitation distribution
qqnorm(CanadianWeather$dailyAv[,,"Precipitation.mm"], datax=TRUE)
# Strongly lognormal?
quantile(CanadianWeather$dailyAv[,,"Precipitation.mm"])
# Can't take logarithms directly,
# because some observations are 0
sum(CanadianWeather$precav==0)
# 27 of 365
# Per help("CanadianWeather"),
sort(unique(diff(sort(unique(CanadianWeather$
                             dailyAv[,,"Precipitation.mm"])))))
# Some repeated numbers, indicating round off problems
sort(unique(diff(sort(round(unique(
     CanadianWeather$dailyAv[,,"Precipitation.mm"]), 7)))))
sort(unique(round(diff(sort(round(unique(
     CanadianWeather$dailyAv[,,"Precipitation.mm"]), 5) )),5)) )
# Obviously, the 0's are problems ... but ignore them
# The real point is that these numbers are all
# to the nearest tenth of a milimeter,
# & the 0 differences are created by anomolies in 'unique'
table(CanadianWeather$dailyAv[,,"Precipitation.mm"])

# help("CanadianWeather") says the 0's were replaced by 0.05 mm
# before computing logarithms
qqnorm(CanadianWeather$dailyAv[,,"log10precip"], datax=TRUE)
# Plausibly a weak annual cycle
# relative to substantial normal noise
# on a log scale

# Conclusion:  Prefer analysis on the log scale
# Back transform to get answers in 'mm' or approx. percent
# (recalling log(1+x) = x if x is small)
dayprecfd <- with(CanadianWeather, smooth.basis(day.5,
     dailyAv[,,"log10precip"], daybasis65,
     fdnames=list("Day", "Station", "log10(mm)"))$fd )
plot(dayprecfd, axes=FALSE)
axisIntervals(1)
axis(2)

# Or create a functional data object with Temp and log10precip together:
CanadianTempPrecip.fd <- with(CanadianWeather, smooth.basis(day.5,
         dailyAv[,,-2], daybasis65)$fd )
str(CanadianTempPrecip.fd)

#  set up plotting arrangements for one and two panel displays allowing
#  for larger fonts

#  Plot temperature curves and values

# This plot would be too busy if we superimposed
# all 35 stations on the same page.

# Therefore, use "index" to make 5 separate plots
# Another alternative would be the 'ask' argument

# Returns invisibly the mean square deviations
# between observed and fdobj=daytempfd
(CanadianWea.MSE <- with(CanadianWeather, plotfit.fd(
   y=dailyAv[,,"Temperature.C"], argvals=day.5,
   fdobj=CanadianTempPrecip.fd[,1], index=1:7, axes=FALSE) ))
axisIntervals(1)
axis(2)

# Same as from a univariate functional data object
(CanadianWea.MSE <- with(CanadianWeather, plotfit.fd(
   y=dailyAv[,,"Temperature.C"], argvals=day.5,
   fdobj=daytempfd, index=1:7, axes=FALSE) ))
axisIntervals(1)
axis(2)

with(CanadianWeather, plotfit.fd(y=dailyAv[,,"Temperature.C"],
   argvals=day.5, fdobj=daytempfd, index=8:14, axes=FALSE) )
axisIntervals(1)
axis(2)

with(CanadianWeather, plotfit.fd(y=dailyAv[,,"Temperature.C"],
   argvals=day.5, fdobj=daytempfd, index=15:21, axes=FALSE) )
axisIntervals(1)
axis(2)

with(CanadianWeather, plotfit.fd(y=dailyAv[,,"Temperature.C"],
   argvals=day.5, fdobj=daytempfd, index=22:28, axes=FALSE) )
axisIntervals(1)
axis(2)

with(CanadianWeather, plotfit.fd(y=dailyAv[,,"Temperature.C"],
   argvals=day.5, fdobj=daytempfd, index=29:35, axes=FALSE) )
axisIntervals(1)
axis(2)

# The smoothing is probably not quite adequate,
# but it's not bad either.

#  Plot residuals for three best fits and three worst fits

#casenames <- CanadianWeather$place
#varnames  <- "Temperature"
#rng       <- dayrange
#index     <- c(1,2,3,33,34,35)
#residual  <- TRUE
#sortwrd   <- TRUE

with(CanadianWeather, plotfit.fd(y=dailyAv[,,"Temperature.C"],
       argvals=day.5, fdobj=daytempfd, index=c(1:3, 33:35),
           nfine=366, residual=TRUE, sortwrd=TRUE, axes=FALSE) )
axisIntervals(1)
axis(2)
# To see the lines individually, plot them 1, 2 or 3 at a time:
#with(CanadianWeather, plotfit.fd(y=dailyAv[,,"Temperature.C",
#   argvals=day.5, fdobj=daytempfd, index=c(1:3, 33:35),
#   nfine=366, residual=TRUE, sortwrd=TRUE, ask=TRUE, col=1, lty=1) )
#with(CanadianWeather, plotfit.fd(y=dailyAV[,,"Temperature.C",
#   argvals=day.5, fdobj=daytempfd, index=c(1:3, 33:35),
#   nfine=366, residual=TRUE, sortwrd=TRUE, ask=TRUE, col=1:3, lty=1) )
# NOTE:  "col=1, lty=1" allows a singly plot per page.

#  Plot precipitation curves and values

CanadianWeather$place
with(CanadianWeather, plotfit.fd(y=dailyAv[,,"log10precip"],
     argvals=day.5, fdobj=dayprecfd, index=1, titles=place,
     axes=FALSE) )
axisIntervals(1)
axis(2)
with(CanadianWeather, plotfit.fd(y=dailyAv[,,"log10precip"],
     argvals=day.5, fdobj=dayprecfd, index=35, titles=place,
     axes=FALSE) )
axisIntervals(1)
axis(2)

#  Assessment: the functions are definitely too rough with
#  this many basis functions, and especially for precip. which
#  has a much higher noise level.

#  These smoothing parameter values probably undersmooth the data,
#  but we can impose further smoothness on the results of analyses

daytempfdSm <- smooth.fdPar(daytempfd, harmaccelLfd365, lambda=10)
dayprecfdSm <- smooth.fdPar(dayprecfd, harmaccelLfd365, lambda=1e5)

# Or
daytempSm <- smooth.fdPar(daytempfdSm, harmaccelLfd365, lambda=10)
str(daytempfdSm)
str(daytempSm)
# fdnames lost in bivariate ...
str(daytempfd)

#  Use function 'plotfit.fd' to plot the temperature data, fit and residuals
#    sorted by size of root mean square error

#  Plot temperature curves and values

with(CanadianWeather, plotfit.fd(y=dailyAv[,,"Temperature.C"],
     argvals=day.5, fdobj=daytempSm, titles=place) )

#  Plot residuals for three best fits and three worst fits

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"],
     day.5, daytempSm, index=c(1,2,3,33,34,35), sortwrd=TRUE,
           residual=TRUE, titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)

#  Plot curves and values only for January
# using rng=c(0, 31)
with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"],
     day.5, daytempSm, rng=c(0,31), titles=place) )


#  plot each pair of functions along with raw data

#par(mfrow=c(1,2), pty="s")
par(mfrow=c(2,1))
for (i in 1:length(CanadianWeather$place) ) {
  with(CanadianWeather, plot(day.5, dailyAv[,i,"Temperature.C"],
       type="p", xlim=c(0, 365), col=1, xlab="Day", ylab="",
         main=paste(place[i],"temperature")) )
  lines.fd(daytempSm[i])
  with(CanadianWeather, plot(day.5, dailyAv[,i,"log10precip"],
       type="p", xlim=c(0, 365), col=1, xlab="Day", ylab="",
         main=paste(place[i],"precipitation")) )
  lines.fd(dayprecfdSm[i])
# Uncomment the following line 'par(ask=TRUE)'
#  to plot the cuves one at a time
# Otherwise, they fly by faster than the eye can see.
   par(ask=TRUE)
}
par(mfrow=c(1,1), pty="m", ask=FALSE)

#  plot all the functions

#par(mfrow=c(1,2), pty="s")
op <- par(mfrow=c(2,1))
plot(daytempfdSm, main="Temperature", axes=FALSE)
axisIntervals(1)
axis(2)
plot(dayprecfdSm, main="Precipitation", axes=FALSE)
axisIntervals(1)
axis(2)
par(op)

#  -------------------------------------------------------------
#                 Choose level of smoothing using
#          the generalized cross-validation criterion
#              with smoothing function smooth.basisPar.
#  -------------------------------------------------------------

wtvec <- rep(1,365)

# set up a saturated basis capable of interpolating the data

daybasis365 <- create.fourier.basis(c(0, 365), 365)

#  --------------------  smooth temperature  ------------------

#  set up range of smoothing parameters in log_10 units

Temp.loglam <- (-5):9
names(Temp.loglam) <- Temp.loglam
Temp.smoothSt <- sapply(Temp.loglam, function(x){
  lam <- 10^x
  smoothList <- with(CanadianWeather, smooth.basisPar(
      argvals=day.5, y=dailyAv[,,"Temperature.C"], fdobj=daybasis365,
      Lfdobj=harmaccelLfd365, lambda=lam) )
  cat(x, "")
  with(smoothList, return(c(loglam=x, df=df, gcv=sum(gcv))))
} )

(Temp.smoothStats <- as.data.frame(t(Temp.smoothSt)))

#library(lattice)
#xyplot(gcv+df~loglam, data=Temp.smoothStats)
#xyplot(gcv+df~loglam, data=Temp.smoothStats,
#       scales=list(y=list(relation="free")))

#op <- par(mfrow=c(1,2), pty="s")
op <- par(mfrow=c(2,1))
with(Temp.smoothStats, {
  plot(loglam, gcv, type="b", cex=1,
     xlab="Log_10 lambda", ylab="GCV Criterion",
     main="Temperature Smoothing")
  plot(loglam, df, type="b",  cex=1,
     xlab="Log_10 lambda", ylab="Degrees of freedom",
     main="Temperature Smoothing", log="y") } )
par(op)

#  Do final smooth with minimum GCV value
# but try other levels also

#lambda   <- 0.01  #  minimum GCV estimate, corresponding to 255 df
#fdParobj <- fdPar(daybasis365, harmaccelLfd, lambda)

#smoothlist. <- with(CanadianWeather, smooth.basis(dayOfYear-0.5,
#                       tempav, fdParobj)$fd )
TempSmooth.01 <- with(CanadianWeather, smooth.basisPar(
     argvals=day.5, y=dailyAv[,, "Temperature.C"], fdobj=daybasis365,
     Lfdobj=harmaccelLfd365, lambda=0.01) )
TempSmooth.1 <- with(CanadianWeather, smooth.basisPar(
     argvals=day.5, y=dailyAv[,, "Temperature.C"], fdobj=daybasis365,
     Lfdobj=harmaccelLfd365, lambda=0.1) )
TempSmooth1 <- with(CanadianWeather, smooth.basisPar(
     argvals=day.5, y=dailyAv[,, "Temperature.C"], fdobj=daybasis365,
     Lfdobj=harmaccelLfd365, lambda=1) )
TempSmooth10 <- with(CanadianWeather, smooth.basisPar(
     argvals=day.5, y=dailyAv[,, "Temperature.C"], fdobj=daybasis365,
     Lfdobj=harmaccelLfd365, lambda=10) )
TempSmooth100 <- with(CanadianWeather, smooth.basisPar(
     argvals=day.5, y=dailyAv[,, "Temperature.C"], fdobj=daybasis365,
     Lfdobj=harmaccelLfd365, lambda=100) )
TempSmooth1000 <- with(CanadianWeather, smooth.basisPar(
     argvals=day.5, y=dailyAv[,, "Temperature.C"], fdobj=daybasis365,
     Lfdobj=harmaccelLfd365, lambda=1000) )
TempSmooth4 <- with(CanadianWeather, smooth.basisPar(
     argvals=day.5, y=dailyAv[,, "Temperature.C"], fdobj=daybasis365,
     Lfdobj=harmaccelLfd365, lambda=1e4) )
TempSmooth5 <- with(CanadianWeather, smooth.basisPar(
     argvals=day.5, y=dailyAv[,, "Temperature.C"], fdobj=daybasis365,
     Lfdobj=harmaccelLfd365, lambda=1e5) )
TempSmooth6 <- with(CanadianWeather, smooth.basisPar(
     argvals=day.5, y=dailyAv[,, "Temperature.C"], fdobj=daybasis365,
     Lfdobj=harmaccelLfd365, lambda=1e6) )
TempSmooth7 <- with(CanadianWeather, smooth.basisPar(
     argvals=day.5, y=dailyAv[,, "Temperature.C"], fdobj=daybasis365,
     Lfdobj=harmaccelLfd365, lambda=1e7) )
TempSmooth8 <- with(CanadianWeather, smooth.basisPar(
     argvals=day.5, y=dailyAv[,, "Temperature.C"], fdobj=daybasis365,
     Lfdobj=harmaccelLfd365, lambda=1e8) )
TempSmooth9 <- with(CanadianWeather, smooth.basisPar(
     argvals=day.5, y=dailyAv[,, "Temperature.C"], fdobj=daybasis365,
     Lfdobj=harmaccelLfd365, lambda=1e9) )

(stderr <- with(TempSmooth.01, sqrt(SSE/(35*(365-df)))))
# 0.26 deg C

#  plot data and fit

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"],
          day.5, TempSmooth.01$fd, titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)

op <- par(xpd=NA, bty="n")
# trim lines at 'device region' not 'plot region'
with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
          TempSmooth.01$fd, index=c(1,35), titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)
lines(TempSmooth.01$fd, lty=1, lwd=2)
title("Canadian Annual Temperature Cycle;  lambda = 0.01")
par(mfrow=c(1,1))

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
          TempSmooth.1$fd, index=c(1,35), titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
          TempSmooth1$fd, index=c(1,35), titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)
lines(TempSmooth1$fd, lty=1, lwd=2)
title("Canadian Annual Temperature Cycle;  lambda = 1")

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
          TempSmooth10$fd, index=c(1,35), titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)
lines(TempSmooth10$fd, lty=1, lwd=2)
title("Canadian Annual Temperature Cycle;  lambda = 10")

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
      TempSmooth100$fd, index=c(1,35), titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)
lines(TempSmooth100$fd, lty=1, lwd=2)
title("Canadian Annual Temperature Cycle;  lambda = 100")

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
       TempSmooth1000$fd, index=c(1,35), titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)
lines(TempSmooth1000$fd, lty=1, lwd=2)
title("Canadian Annual Temperature Cycle;  lambda = 1000")

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
      TempSmooth4$fd, index=c(1,35), titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)
lines(TempSmooth4$fd, lty=1, lwd=2)
title("Canadian Annual Temperature Cycle;  lambda = 1e4")

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
          TempSmooth5$fd, index=c(1,35), titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)
lines(TempSmooth5$fd, lty=1, lwd=2)
title("Canadian Annual Temperature Cycle;  lambda = 1e5")

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
      TempSmooth6$fd, index=c(1,35), titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)
lines(TempSmooth6$fd, lty=1, lwd=2)
title("Canadian Annual Temperature Cycle;  lambda = 1e6")
# Smoothing with lambda = 1e6 follows the main pattern
# but NOT the fine detail.
# Q:  Is that fine detail real or serial dependence
# IGNORED by gcv?  If the latter, we should use
# lambda = 1e6 over the 'gcv' min of 0.01.

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
      TempSmooth7$fd, index=c(1,35), titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)
lines(TempSmooth7$fd, lty=1, lwd=2)
title("Canadian Annual Temperature Cycle;  lambda = 1e7")
# possibly acceptable

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
      TempSmooth8$fd, index=c(1,35), titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)
lines(TempSmooth8$fd, lty=1, lwd=2)
title("Canadian Annual Temperature Cycle;  lambda = 1e8")
# subtle but clear oversmoothing

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
       TempSmooth9$fd, index=c(1,35), titles=place, axes=FALSE) )
axisIntervals(1)
axis(2)
lines(TempSmooth9$fd, lty=1, lwd=2)
title("Canadian Annual Temperature Cycle;  lambda = 1e9")
# lambda=1e9 is oversmoothing, blatently obvious

lvl.1 <- (-10:10)/10
contour(cor.fd(weeks, TempSmooth1000$fd), levels=lvl.1)
# NOT as smooth as Fig. 2.4, FDA
op <- par(mfrow=c(2,2))
contour(weeks, weeks, cor.fd(weeks, TempSmooth4$fd), levels=lvl.1,
        xlab="Average daily Temperature (C)",
        ylab="Average daily Temperature (C)",
        axes=FALSE, main="Temperature correlations, smoothing = 1e4")
axisIntervals(1, labels=monthLetters)
axisIntervals(2, labels=monthLetters)
contour(weeks, weeks, cor.fd(weeks, TempSmooth5$fd), levels=lvl.1,
        xlab="Average daily Temperature (C)",
        ylab="Average daily Temperature (C)",
        axes=FALSE, main="Temperature correlations, smoothing = 1e5")
axisIntervals(1, labels=monthLetters)
axisIntervals(2, labels=monthLetters)
contour(weeks, weeks, cor.fd(weeks, TempSmooth6$fd), levels=lvl.1,
        xlab="Average daily Temperature (C)",
        ylab="Average daily Temperature (C)",
        axes=FALSE, main="Temperature correlations, smoothing = 1e6")
axisIntervals(1, labels=monthLetters)
axisIntervals(2, labels=monthLetters)
contour(weeks, weeks, cor.fd(weeks, TempSmooth7$fd), levels=lvl.1,
        xlab="Average daily Temperature (C)",
        ylab="Average daily Temperature (C)",
        axes=FALSE, main="Temperature correlations, smoothing = 1e7")
axisIntervals(1, labels=monthLetters)
axisIntervals(2, labels=monthLetters)
par(op)
# Figure 2.4 looks closest to lambda = 1e6.

#  --------------------  smooth precipitation  ------------------

Prec.loglam <- -5:9
names(Prec.loglam) <- Prec.loglam
Prec.smoothSt <- sapply(Prec.loglam, function(x){
  lam <- 10^x
  smoothList <- with(CanadianWeather, smooth.basisPar(
      argvals=day.5, y=dailyAv[,,"log10precip"], fdobj=daybasis365,
      Lfdobj=harmaccelLfd365, lambda=lam) )
  cat(x, "")
  with(smoothList, return(c(loglam=x, df=df, gcv=sum(gcv))))
} )

(Prec.smoothStats <- as.data.frame(t(Prec.smoothSt)))

op <- par(mfrow=c(2,1))
with(Prec.smoothStats[8:16,], {
  plot(loglam, gcv, type="b", cex=1, log="y",
     xlab="Log_10 lambda", ylab="GCV Criterion",
     main="Precipitation Smoothing")
  plot(loglam, df, type="b",  cex=1,
     xlab="Log_10 lambda", ylab="Degrees of freedom",
     main="Precipitation Smoothing", log="y") } )
par(op)

# lambda = 1e6 minimizes gcv, df = 12.3

#  Do final smooth with minimum GCV value

# Previous note:
#lambda   <- 1e7  #  minimum GCV estimate, corresponding to 255 df
# The df can NOT be correct with this lambda

PrecSmooth6 <- with(CanadianWeather, smooth.basisPar(
    argvals=day.5, y=dailyAv[,,"log10precip"],
    fdobj=daybasis365, Lfdobj=harmaccelLfd365, lambda=1e6) )

(stderr <- with(PrecSmooth6, sqrt(SSE/(35*(365-df)))))
# 0.198 vs. previous annotation of 0.94 ???

class(PrecSmooth6)
sapply(PrecSmooth6, class)
# Looks oversmoothed relative to Figure 2.4, FDA

PrecSmooth0 <- with(CanadianWeather, smooth.basisPar(
    argvals=day.5, y=dailyAv[,,"log10precip"],
    fdobj=daybasis365, Lfdobj=harmaccelLfd365, lambda=1) )
contour(cor.fd(weeks, PrecSmooth0$fd), levels=lvl.1)
# Way undersmoothed

PrecSmooth3 <- with(CanadianWeather, smooth.basisPar(
    argvals=day.5, y=dailyAv[,,"log10precip"],
    fdobj=daybasis365, Lfdobj=harmaccelLfd365, lambda=1e3) )
# Still undersmoothed but not as bad as lambda=1

PrecSmooth4 <- with(CanadianWeather, smooth.basisPar(
    argvals=day.5, y=dailyAv[,,"log10precip"],
    fdobj=daybasis365, Lfdobj=harmaccelLfd365, lambda=1e4) )

PrecSmooth5 <- with(CanadianWeather, smooth.basisPar(
    argvals=day.5, y=dailyAv[,,"log10precip"],
    fdobj=daybasis365, Lfdobj=harmaccelLfd365, lambda=1e5) )



op <- par(mfrow=c(2,2))
contour(weeks, weeks, cor.fd(weeks, PrecSmooth3$fd), levels=lvl.1,
        xlab="Average daily precipitation (mm)",
        ylab="Average daily precipitation (mm)",
        axes=FALSE, main="Precipitation correlations, smoothing 1e3")
axisIntervals(1, labels=monthLetters)
axisIntervals(2, labels=monthLetters)

contour(weeks, weeks, cor.fd(weeks, PrecSmooth4$fd), levels=lvl.1,
        xlab="Average daily precipitation (mm)",
        ylab="Average daily precipitation (mm)",
        axes=FALSE, main="Precipitation correlations, smoothing 1e4")
axisIntervals(1, labels=monthLetters)
axisIntervals(2, labels=monthLetters)

contour(weeks, weeks, cor.fd(weeks, PrecSmooth5$fd), levels=lvl.1,
        xlab="Average daily precipitation (mm)",
        ylab="Average daily precipitation (mm)",
        axes=FALSE, main="Precipitation correlations, smoothing 1e5")
axisIntervals(1, labels=monthLetters)
axisIntervals(2, labels=monthLetters)

contour(weeks, weeks, cor.fd(weeks, PrecSmooth6$fd), levels=lvl.1,
        xlab="Average daily precipitation (mm)",
        ylab="Average daily precipitation (mm)",
        axes=FALSE, main="Precipitation correlations, smoothing 1e6")
axisIntervals(1, labels=monthLetters)
axisIntervals(2, labels=monthLetters)

par(op)
# lambda = 1e6 looks the closest to Fig. 2.4, FDA

# cross correlations
contour(weeks, weeks,
   cor.fd(weeks, TempSmooth6$fd, weeks, PrecSmooth6$fd), levels=lvl.1,
        xlab="Average daily Temperature (C)",
        ylab="Average daily precipitation (mm)",
        axes=FALSE, main="Canadian Weather Correlations")
axisIntervals(1)
axisIntervals(2)

# CONCLUSIONS FROM CROSS CORRELATIONS:
# Places where January is warmer are quite likely to have more rain (r=.8).
# Places where June is warmer are only moderately to have more rain (r=.4).


#  plot data and fit

#par(mfrow=c(1,1), pty="m")
with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"],
          day.5, TempSmooth6$fd, titles=place) )
with(CanadianWeather, plotfit.fd(dailyAv[,,"log10precip"],
          day.5, PrecSmooth6$fd, titles=place) )
with(CanadianWeather, plotfit.fd(dailyAv[,,"log10precip"],
          day.5, PrecSmooth6$fd, titles=place, index=1:2) )
with(CanadianWeather, plotfit.fd(dailyAv[,,"log10precip"],
          day.5, PrecSmooth6$fd, titles=place, index=c(2, 35)) )

#with(CanadianWeather, plotfit.fd(dailyAv[,,"log10precip"], day.5,
#          PrecSmooth6$fd, titles=place, lty=1, col=1, ask=TRUE) )

#  Assessment: the temperature curves are still pretty rough,
#  although the data themselves show that there are very
#  high frequency effects in the mean temperature, especially
#  early in the year.
# From the indivual plots (ask=TRUE, lty=1, col=1),
# Ramsay & Silverman felt that the precip. curves may be
# oversmoothed for some weather stations (?)

#  smooth precipitation in Prince Rupert

#PRprecfd <- smooth.basis(daytime, CanadianWeather$precav[,29], fdParobj)$fd
PRprecfd <- smooth.basisPar(day.5, CanadianWeather$dailyAv[,29,"log10precip"],
                            PrecSmooth6$fd, harmaccelLfd365, lambda=1e6)

#PRprecvec <- eval.fd(day.5, PRprecfd)
PRprecvec <- predict(PRprecfd, day.5)

plot(day.5, CanadianWeather$dailyAv[,29,"log10precip"], type="p",
     xlab="Day", ylab="Precipitation (mm)", main="Prince Rupert")
lines(day.5, PRprecvec, lwd=2)

#  ----------------------------------------------------------------------
#                Descriptive Statistics Functions
#  ----------------------------------------------------------------------

#  ---------  create fd objects for temp. and prec. ---------------

# First check the distribution
qqnorm(CanadianWeather$dailyAv[,,"Temperature.C"], datax=TRUE)
# Short tailed distribution = strong annual cycle +
# plausibly modest normal error

# Recreate daytempfd, dayprecfd created above
daytempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],
         daybasis65, fdnames=list("Day", "Station", "Deg C"))$fd

dayprecfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"log10precip"],
         daybasis65, fdnames=list("Day", "Station", "log10(mm)"))$fd

#  --  compute and plot mean and standard deviation of temperature -------

(tempmeanfd  <- mean.fd(daytempfd))
(tempstdvfd  <- sd.fd(daytempfd))

#op <- par(mfrow=c(1,2), pty="s")
#op <- par(mfrow=c(2,1))
#plot(tempmeanfd,               main="Mean")
#plot(tempstdvfd, ylim=c(0,10), main="Standard Deviation")
#par(op)
op <- par(mfrow=c(2,1))
plot(tempmeanfd,               main="Mean")
plot(tempstdvfd, main="Standard Deviation", log="y")
par(op)

#  --  plot the temperature variance-covariance bivariate function  ----

str(tempvarbifd <- var.fd(daytempfd))
str(tempvarmat  <- eval.bifd(weeks,weeks,tempvarbifd))
# dim(tempvarmat)= c(53, 53)

op <- par(mfrow=c(1,2), pty="s")
#contour(tempvarmat, xlab="Days", ylab="Days")
contour(weeks, weeks, tempvarmat,
        xlab="Temperature by day",
        ylab="Temperature by day",
        main=paste("Variance function across locations\n",
          "for Canadian Anual Temperature Cycle"),
        cex.main=0.8, axes=FALSE)
axisIntervals(1, atTick1=seq(0, 365, length=5), atTick2=NA,
            atLabels=seq(1/8, 1, 1/4)*365,
            labels=paste("Q", 1:4) )
axisIntervals(2, atTick1=seq(0, 365, length=5), atTick2=NA,
            atLabels=seq(1/8, 1, 1/4)*365,
            labels=paste("Q", 1:4) )
#persp(tempvarmat,xlab="Days", ylab="Days", zlab="Covariance")
persp(weeks, weeks, tempvarmat,
      xlab="Days", ylab="Days", zlab="Covariance")
mtext("Temperature Covariance", line=-4, outer=TRUE)
par(op)

#  --  plot the temperature correlation function  ---

tempstddev <- sqrt(diag(tempvarmat))
tempcormat <- tempvarmat/outer(tempstddev,tempstddev)
tempcormat <- tempvarmat/(tempstddev %o% tempstddev)

op <- par(mfrow=c(1,2), pty="s")
contour(weeks, weeks, tempcormat, xlab="Days", ylab="Days", axes=FALSE)
axisIntervals(1, atTick1=seq(0, 365, length=5), atTick2=NA,
            atLabels=seq(1/8, 1, 1/4)*365,
            labels=paste("Q", 1:4) )
axisIntervals(2, atTick1=seq(0, 365, length=5), atTick2=NA,
            atLabels=seq(1/8, 1, 1/4)*365,
            labels=paste("Q", 1:4) )
persp(weeks, weeks, tempcormat,
      xlab="Days", ylab="Days", zlab="Covariance")
mtext("Temperature Correlation", line=-4, outer=TRUE)
par(op)
#  --  plot the precipitation variance-covariance bivariate function  ----

#weeks <- seq(0,365,length=53)
precvarbifd <- var.fd(dayprecfd)
precvarmat  <- eval.bifd(weeks,weeks,precvarbifd)

op <- par(mfrow=c(2,2), pty="s")
contour(weeks, weeks, precvarmat, xlab="Days", ylab="Days", axes=FALSE)
axisIntervals(1, atTick1=seq(0, 365, length=5), atTick2=NA,
            atLabels=seq(1/8, 1, 1/4)*365,
            labels=paste("Q", 1:4) )
axisIntervals(2, atTick1=seq(0, 365, length=5), atTick2=NA,
            atLabels=seq(1/8, 1, 1/4)*365,
            labels=paste("Q", 1:4) )

persp(weeks, weeks, precvarmat,
      xlab="Days", ylab="Days", zlab="Covariance")
mtext("Precipitation Covariance", line=-4, outer=TRUE)

#  --  plot the precipitation correlation function  ---

precstddev <- sqrt(diag(precvarmat))
preccormat <- precvarmat/(precstddev %o% precstddev)

contour(preccormat, xlab="Days", ylab="Days", axes=FALSE)
axisIntervals(1, atTick1=seq(0, 365, length=5), atTick2=NA,
            atLabels=seq(1/8, 1, 1/4)*365,
            labels=paste("Q", 1:4) )
axisIntervals(2, atTick1=seq(0, 365, length=5), atTick2=NA,
            atLabels=seq(1/8, 1, 1/4)*365,
            labels=paste("Q", 1:4) )
mtext("Precipitation Correlation")
persp(preccormat,
      xlab="Days", ylab="Days", zlab="Covariance")
mtext("Precipitation Correlation")
par(op)

#  -----  compute and plot the covariance between temp. and prec.  --

covbifd <- var.fd(daytempfd, dayprecfd)
covmat  <- eval.bifd(weeks,weeks,covbifd)

contour(weeks, weeks, covmat, xlab="Weeks", ylab="Weeks", axes=FALSE)
axisIntervals(1)
axisIntervals(2)

persp(weeks, weeks, covmat, xlab="Weeks", ylab="Weeks", zlab="Covariance")
mtext("Temperature-Precipitation Covariance", line=-4, outer=TRUE)

#  -----  compute and plot the correlation between temp. and prec.  --

cormat  <- covmat/(tempstddev %o% precstddev)

contour(weeks, weeks, cormat, xlab="Weeks", ylab="Weeks", axes=FALSE)
axisIntervals(1)
axisIntervals(2)

persp(cormat, xlab="Weeks", ylab="Weeks", zlab="Correlation")
mtext("Temperature-Precipitation Correlation", line=-4, outer=TRUE)


#  -----------------------------------------------------------------------
#               PCA of temperatures with varimax rotation
#  -----------------------------------------------------------------------

harmfdPar     <- fdPar(daybasis65, harmaccelLfd365, 1e5)

daytemppcaobj <- pca.fd(daytempfd, nharm=4, harmfdPar)

daytemppcaobjVM <- varmx.pca.fd(daytemppcaobj)
str(daytemppcaobjVM)
dimnames(daytemppcaobjVM$scores)[[2]] <- paste("PCA", 1:4, sep=".")
round(daytemppcaobjVM$scores)

#  plot harmonics

par(mfrow=c(1,1), pty="m")
plot.pca.fd(daytemppcaobjVM)

op <- par(mfrow=c(2,2), pty="m")
plot.pca.fd(daytemppcaobjVM)
par(op)

#  plot log eigenvalues

daytempeigvals <- daytemppcaobjVM[[2]]

plot(1:20, log10(daytempeigvals[1:20]), type="b",
     xlab="Eigenvalue Number", ylab="Log 10 Eigenvalue")
abline(lsfit(5:20, log10(daytempeigvals[5:20])), lty=2)

#  plot factor scores

harmscr <- daytemppcaobjVM[[3]]

plot(harmscr[,1], harmscr[,2], xlab="Harmonic 1", ylab="Harmonic 2")
text(harmscr[,1], harmscr[,2], CanadianWeather$place, col=4)

plot(harmscr[,3], harmscr[,4], xlab="Harmonic 3", ylab="Harmonic 4")
text(harmscr[,3], harmscr[,4], CanadianWeather$place, col=4)

#  ------------------------------------------------------------------
#               Functional linear models
#  ------------------------------------------------------------------

#  ---------------------------------------------------------------
#             Predicting temperature from climate zone
#  ---------------------------------------------------------------

#  return data to original ordering

# CanadianWeather <- vector("list", 0)
# CanadianWeather$tempav <- matrix(scan("../data/dailtemp.txt",0), 365, 35)
# CanadianWeather$precav <- matrix(scan("../data/dailprec.txt",0), 365, 35)

#  set up a smaller basis using only 65 Fourier basis functions
#  to save some computation time

#smallnbasis <- 65
#smallbasis  <- create.fourier.basis(c(0, 365), smallnbasis)
smallbasis  <- create.fourier.basis(c(0, 365), 65)

tempfd      <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],
                       smallbasis)$fd

smallbasismat <- eval.basis(day.5, smallbasis)
y2cMap <- solve(crossprod(smallbasismat), t(smallbasismat))

#  names for climate zones

zonenames <- c("Canada  ",
               "Atlantic", "Pacific ", "Contintal", "Arctic  ")

#  indices for (weather stations in each of four climate zones

index = 1:35


atlindex <- index[CanadianWeather$region == "Atlantic"]
pacindex <- index[CanadianWeather$region == "Pacific"]
conindex <- index[CanadianWeather$region == "Continental"]
artindex <- index[CanadianWeather$region == "Arctic"]

#  Set up a design matrix having a column for (the grand mean, and
#    a column for (each climate zone effect. Add a dummy contraint
#    observation

zmat <- matrix(0,35,5)
zmat[        ,1] <- 1
zmat[atlindex,2] <- 1
zmat[pacindex,3] <- 1
zmat[conindex,4] <- 1
zmat[artindex,5] <- 1

#  attach a row of 0, 1, 1, 1, 1 to force zone
#  effects to sum to zero, and define first regression
#  function as grand mean for (all stations

z36    <- matrix(1,1,5)
z36[1] <- 0
zmat   <- rbind(zmat, z36)

#  revise YFDOBJ by adding a zero function

coef   <- tempfd$coefs
str(coef)
# add a 0 column # 36 to coef
coef36 <- cbind(coef,matrix(0,65,1))
tempfd$coefs <- coef36

p <- 5
xfdlist <- vector("list",p)
for (j in 1:p) xfdlist[[j]] <- zmat[,j]

#  set up the basis for (the regression functions

nbetabasis <- 11
betabasis  <- create.fourier.basis(c(0, 365), nbetabasis)

#  set up the functional parameter object for (the regression fns.

betafd    <- fd(matrix(0,nbetabasis,1), betabasis)
estimate  <- TRUE
lambda    <- 0
betafdPar <- fdPar(betafd, harmaccelLfd365, lambda, estimate)

betalist <- vector("list",p)
for (j in 1:p) betalist[[j]] <- betafdPar

#  compute regression coefficient functions and
#  predicted functions

fRegressList <- fRegress(tempfd, xfdlist, betalist)


#  plot regression functions

betaestlist <- fRegressList$betaestlist
par(mfrow=c(3,2))
for (j in 1:p) {
	betaestParfdj <- betaestlist[[j]]
	plot(betaestParfdj$fd, xlab="Day", ylab="Temp.",
	     main=zonenames[j])
}

#  plot predicted functions

yhatfdobj <- fRegressList$yhatfdobj
plot(yhatfdobj,main='Predicted Temperature',)

#  compute residual matrix and get covariance of residuals

#yhatmat  <- eval.fd(day.5, yhatfdobj)
yhatmat  <- predict(yhatfdobj, day.5)
ymat     <- eval.fd(day.5, tempfd)
temprmat <- ymat[,1:35] - yhatmat[,1:35]
SigmaE   <- var(t(temprmat))

#  plot covariance surface for errors

par(mfrow=c(1,1))
contour(SigmaE, xlab="Day", ylab="Day")
lines(c(0, 365), c(0, 365),lty=4)

#  plot standard deviation of errors

par(mfrow=c(1,1), mar=c(5,5,3,2), pty="m")
stddevE <- sqrt(diag(SigmaE))
plot(day.5, stddevE, type="l",
     xlab="Day", ylab="Standard error (deg C)")

#  Repeat regression, this time outputting results for
#  confidence intervals

stderrList <- fRegress.stderr(fRegressList, y2cMap, SigmaE)

betastderrlist <- stderrList$betastderrlist

#  plot regression function standard errors

op <- par(mfrow=c(2,3), pty="s")
for (j in 1:p) {
	betastderrj <- eval.fd(day.5, betastderrlist[[j]])
	plot(day.5, betastderrj,
	        type="l",lty=1, xlab="Day", ylab="Reg. Coeff.",
	        main=zonenames[j])
}
par(op)

#  plot regression functions with confidence limits

op <- par(mfrow=c(3,2))
for (j in 1:p) {
	betafdParj  <- betaestlist[[j]]
	betafdj     <- betafdParj$fd
	betaj       <- eval.fd(day.5, betafdj)
	betastderrj <- eval.fd(day.5, betastderrlist[[j]])
	matplot(day.5, cbind(betaj, betaj+2*betastderrj, betaj-2*betastderrj),
	        type="l",lty=c(1,4,4), xlab="Day", ylab="Reg. Coeff.",
	        main=zonenames[j])
}
par(op)



# Now a couple of permutation tests

# permutation t-test between atlantic and pacific

t.res = tperm.fd(tempfd[atlindex],tempfd[pacindex])


# instead, we'll try a permutation F-test for the regression

F.res = Fperm.fd(tempfd, xfdlist, betalist,cex.axis=1.5,cex.lab=1.5)




#  -----------------------------------------------------------------------
#         predict log precipitation from climate zone and temperature
#  -----------------------------------------------------------------------

#  Be sure to run previous analysis predicting temperature from
#  climate zone before running this example.

#  set up functional data object for log precipitation

precfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Precipitation.mm"],
                  smallbasis)$fd

logprecmat <- log10(eval.fd(day.5, precfd))

lnprecfd <- smooth.basis(day.5, logprecmat, smallbasis)$fd
lnprecfd$fdnames[[1]] <- "Days"
lnprecfd$fdnames[[2]] <- "Station"
lnprecfd$fdnames[[3]] <- "log.{10} mm"

#  plot log precipitation functions

par(mfrow=c(1,1), pty="m")
plot(lnprecfd)
title("Log Precipitation Functions")

#  revise LOGPREDFD by adding a zero function

coef   <- lnprecfd$coefs
nbasis <- smallbasis$nbasis
coef36 <- cbind(coef,matrix(0,nbasis,1))
lnprecfd$coefs <- coef36

#  set up the XFDLIST list

p <- 6
xfdlist <- vector("list",p)

#  load first five members with columns of design matrix

for (j in 1:5) xfdlist[[j]] <- zmat[,j]

#  set up a FD object for (temperature residuals

lambda     <- 1e5
fdParobj   <- fdPar(smallbasis, harmaccelLfd365, lambda)
smoothList <- smooth.basis(day.5, temprmat, fdParobj)
temprfdobj <- smoothList$fd

#  plot temperature residuals

par(mfrow=c(1,1), pty="m")
plot(temprfdobj)

#  extend temperature residual functions to include
#  zero function

coef   <- temprfdobj$coefs
nbasis <- dim(coef)[1]
coef36 <- cbind(coef,matrix(0,nbasis,1))
temprfdobj$coefs <- coef36

#  add TEMPRFDOBJ to the set of predictors

xfdlist[[6]]  <- temprfdobj
betalist[[6]] <- betafdPar

#  set up the basis for (the regression functions

nbetabasis <- 13
betabasis  <- create.fourier.basis(c(0, 365), nbetabasis)

#  set up the functional parameter object for (the regression fns.

betafd    <- fd(matrix(0,nbetabasis,p), betabasis)
estimate  <- TRUE
lambda    <- 0
betafdPar <- fdPar(betafd, harmaccelLfd365, lambda, estimate)
for (j in 1:p) betalist[[j]] <- betafdPar

#  compute regression coefficient functions and
#  predicted functions

fRegressList <- fRegress(lnprecfd, xfdlist, betalist)

betaestlist <- fRegressList$betaestlist
yhatfdobj   <- fRegressList$yhatfdobj

#  plot regression functions

prednames <- c(zonenames, "tempres ")
op <- par(mfrow=c(2,3),pty="s")
for (j in 1:p) {
	betaParfdj <- betaestlist[[j]]
	betafdj    <- betaParfdj$fd
    plot(betafdj)
    title(prednames[j])
}
par(op)

#  plot predicted functions

par(mfrow=c(1,1), pty="m")
plot(yhatfdobj)

#  compute residual matrix and get covariance of residuals

#yhatmat    <- eval.fd(day.5, yhatfdobj)
yhatmat    <- predict(yhatfdobj, day.5)
ymat       <- eval.fd(day.5, lnprecfd)
lnprecrmat <- ymat[,1:35] - yhatmat[,1:35]
SigmaE     <- var(t(lnprecrmat))

contour(SigmaE)

#  repeat regression analysis to get confidence intervals

stderrList <- fRegress.stderr(fRegressList, y2cMap, SigmaE)

betastderrlist <- stderrList$betastderrlist

#  plot regression functions

prednames <- c(zonenames, "tempres ")

#  plot regression function standard errors

op <- par(mfrow=c(2,3), pty="s")
for (j in 1:p) {
	betastderrfdj <- betastderrlist[[j]]
	betastderrj <- eval.fd(day.5, betastderrfdj)
	plot(day.5, betastderrj,
	        type="l",lty=1, xlab="Day", ylab="Reg. Coeff.",
	        main=prednames[j])
}
par(op)

#  plot regression functions with confidence limits

op <- par(mfrow=c(2,3), pty="s")
for (j in 1:p) {
	betafdParj  <- betaestlist[[j]]
	betafdj     <- betafdParj$fd
	betaj       <- eval.fd(day.5, betafdj)
	betastderrfdj <- betastderrlist[[j]]
	betastderrj   <- eval.fd(day.5, betastderrfdj)
	matplot(day.5, cbind(betaj, betaj+2*betastderrj, betaj-2*betastderrj),
	        type="l",lty=c(1,4,4), xlab="Day", ylab="Reg. Coeff.",
	        main=prednames[j])
}
par(op)



#  ---------------------------------------------------------------
#      log annual precipitation predicted by temperature profile
#  ---------------------------------------------------------------

#  set up log10 total precipitation

annualprec <- log10(apply(CanadianWeather$dailyAv[,,"Precipitation.mm"],
                          2,sum))

#  set up a smaller basis using only 65 Fourier basis functions
#  to save some computation time

smallnbasis <- 65
smallbasis  <- create.fourier.basis(c(0, 365), smallnbasis)
tempfd      <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],
                       smallbasis)$fd

smallbasismat <- eval.basis(day.5, smallbasis)
y2cMap <- solve(crossprod(smallbasismat)) %*% t(smallbasismat)

#  set up the covariates, the first the constant, and the second
#  temperature

p <- 2
constantfd <- fd(matrix(1,1,35), create.constant.basis(c(0, 365)))

xfdlist <- vector("list",2)
xfdlist[[1]] <- constantfd
xfdlist[[2]] <- tempfd[1:35]

#  set up the functional parameter object for (the regression fns.
#  the smoothing parameter for (the temperature function
#  is obviously too small here, and will be revised below by
#  using cross-validation.

betalist   <- vector("list",2)
#  set up the first regression function as a constant
betabasis1 <- create.constant.basis(c(0, 365))
betafd1    <- fd(0, betabasis1)
betafdPar1 <- fdPar(betafd1)
betalist[[1]] <- betafdPar1
#  set up the second with same basis as for (temperature
#  35 basis functions would permit a perfect fit to the data
nbetabasis  <- 35
betabasis2  <- create.fourier.basis(c(0, 365), nbetabasis)
betafd2     <- fd(matrix(0,nbetabasis,1), betabasis2)
lambda      <- 10
betafdPar2  <- fdPar(betafd2, harmaccelLfd365, lambda)
betalist[[2]] <- betafdPar2

#  carry out the regression analysis

fRegressList <- fRegress(annualprec, xfdlist, betalist)

betaestlist   <- fRegressList$betaestlist

annualprechat <- fRegressList$yhatfdobj

#  constant term

alphafdPar <- betaestlist[[1]]

alphafdPar$fd$coefs

#  plot the coefficient function for (temperature

betafdPar <- betaestlist[[2]]
betafd    <- betafdPar$fd
par(mfrow=c(1,1), pty="m")
plot(betafd)
title("Regression coefficient for temperature")

#  plot the fit

plot (annualprechat, annualprec, type="p")
lines(annualprechat, annualprechat, lty=4)

#  compute cross-validated SSE"s for a range of smoothing parameters

loglam <- seq(5,15,0.5)
nlam   <- length(loglam)
SSE.CV <- matrix(0,nlam,1)
for (ilam in 1:nlam) {
    lambda       <- 10^loglam[ilam]
    betalisti    <- betalist
    betafdPar2   <- betalisti[[2]]
    betafdPar2$lambda <- lambda
    betalisti[[2]] <- betafdPar2
    SSE.CV[ilam]   <- fRegress.CV(annualprec, xfdlist, betalisti)
    print(c(ilam, loglam[ilam], SSE.CV[ilam]))
}

plot(loglam, SSE.CV, type="b",
	  xlab="log smoothing parameter lambda",
	  ylab="Cross-validation score",cex.lab=1.5,cex.axis=1.5,
	  lwd=2)

#  analysis with minimum CV smoothing

lambda        <- 10^12.5
betafdPar2    <- fdPar(betafd2, harmaccelLfd365, lambda)
betalist[[2]] <- betafdPar2

#  carry out the regression analysis

fRegressList <- fRegress(annualprec, xfdlist, betalist)

betaestlist   <- fRegressList$betaestlist

annualprechat <- fRegressList$yhatfdobj

#  constant term

alphafdPar <- betaestlist[[1]]

alphafdPar$fd$coefs

#  plot the coefficient function for (temperature

betafdPar <- betaestlist[[2]]
betafd    <- betafdPar$fd

par(mfrow=c(1,1), pty="m")
plot(betafd)
title("Regression coefficient for temperature")

#  plot the fit

par(mfrow=c(1,1), pty="m")
plot (annualprechat, annualprec, type="p",cex.lab=1.5,cex.axis=1.5)
lines(annualprechat, annualprechat, lty=2)

#  compute squared multiple correlation

covmat <- var(cbind(annualprec, annualprechat))
Rsqrd <- covmat[1,2]^2/(covmat[1,1]*covmat[2,2])
#   0.7540

#  compute SigmaE

resid  <- annualprec - annualprechat
SigmaE <- sum(resid^2)/(35-fRegressList$df)
SigmaE <- SigmaE*diag(rep(1,35))

#  recompute the analysis to get confidence limits

stderrList <- fRegress.stderr(fRegressList, NULL, SigmaE)

betastderrlist <- stderrList$betastderrlist

#  constant  coefficient standard error:

betastderr1 <- betastderrlist[[1]]
betastderr1$coefs

#  plot the temperature coefficient function

betafdParj      <- betaestlist[[2]]
betafd          <- betafdParj$fd
betastderrfd    <- betastderrlist[[2]]
betavec         <- eval.fd(day.5, betafd)
betastderrvec   <- eval.fd(day.5, betastderrfd)

betaplotmat <- cbind(betavec, betavec+2*betastderrvec,
                              betavec-2*betastderrvec)

matplot(day.5, betaplotmat, type="l", lty=c(1,4,4),
        xlab="Day", ylab="Temperature Reg. Coeff.",lwd=2,
	 cex.lab=1.5,cex.axis=1.5)
lines(c(0, 365),c(0,0),lty=2)


# We can F-test this as well

F.res2 = Fperm.fd(annualprec, xfdlist, betalist)



#  ---------------------------------------------------------------
#         predict log precipitation from temperature
#  ---------------------------------------------------------------

#  set up a smaller basis using only 65 Fourier basis functions
#  to save some computation time

smallnbasis <- 65
smallbasis  <- create.fourier.basis(c(0, 365), smallnbasis)

tempfd      <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],
                       smallbasis)$fd

#  change 0's to 0.05 mm in precipitation data

prectmp <- CanadianWeather$dailyAv[,,"Precipitation.mm"]
for (j in 1:35) {
    index <- prectmp[,j] == 0
    prectmp[index,j] <- 0.05
}

#  work with log base 10 precipitation

logprecmat <- log10(prectmp)

#  set up functional data object for log precipitation

lnprecfd <- smooth.basis(day.5, logprecmat, smallbasis)$fd
lnprecfdnames <- vector("list",3)
lnprecfdnames[[1]] <- "Days"
lnprecfdnames[[2]] <- "Station"
lnprecfdnames[[3]] <- "log.{10} mm"
lnprecfd$fdnames <- lnprecfdnames

#  plot precipitation functions

par(mfrow=c(1,1), pty="m")
plot(lnprecfd)
title("Log Precipitation Functions")

#  set up smoothing levels for (s (xLfd) and for (t (yLfd)

xLfdobj <- harmaccelLfd365
yLfdobj <- harmaccelLfd365
xlambda <- 1e9
ylambda <- 1e7

#  compute the linear model

wtvec <- matrix(1,35,1)

linmodstr <- linmod(daytempfd, lnprecfd, wtvec,
                    xLfdobj, yLfdobj, xlambda, ylambda)

afd <- linmodstr$alphafd   #  The intercept function
bfd <- linmodstr$regfd     #  The bivariate regression function

#  plot the intercept function

plot(afd, xlab="Day", ylab="Intercept function")

#  plot the regression function as a surface

bfdmat <- eval.bifd(weeks, weeks, bfd)

persp(weeks, weeks, bfdmat, xlab="Day(t)", ylab="Day(s)")

#  Get fitted functions

lnprechatfd <- linmodstr$yhatfd

# Compute mean function as a benchmark for (comparison

lnprecmeanfd <- mean(lnprecfd)
lnprechat0   <- eval.fd(weeks, lnprecmeanfd)

#  Plot actual observed, fitted, and mean log precipitation for
#      each weather station,

lnprecmat    <- eval.fd(weeks, lnprecfd)
lnprechatmat <- eval.fd(weeks, lnprechatfd)
plotrange    <- c(min(lnprecmat),max(lnprecmat))
#par(ask=TRUE)
for (i in 1:35) {
    lnpreci    <- eval.fd(lnprecfd[i],    weeks)
    lnprechati <- eval.fd(lnprechatfd[i], weeks)
    SSE <- sum((lnpreci-lnprechati)^2)
    SSY <- sum((lnpreci-lnprechat0)^2)
    RSQ <- (SSY-SSE)/SSY
    plot(weeks, lnprecmat[,i], type="l", lty=4, ylim=plotrange,
         xlab="Day", ylab="Log Precipitation",
         main=paste(CanadianWeather$place[i],"  R^2 =",round(RSQ,3)))
    lines(weeks, lnprechatmat[,i])
}
par(ask=FALSE)

#  -------------------------------------------------------------------
#              Smooth Vancouver's precipitation with a
#                       positive function.
#  -------------------------------------------------------------------

#  select Vancouver's precipitation

index <- (1:35)[CanadianWeather$place == "Vancouver"]
VanPrec  <- CanadianWeather$dailyAv[,index, "Precipitation.mm"]

#  smooth the data using 65 basis functions

lambda    <- 1e4
fdParobj  <- fdPar(smallbasis, harmaccelLfd365, lambda)
VanPrecfd <- smooth.basis(day.5, VanPrec, fdParobj)$fd

#  Plot temperature curves and values

plotfit.fd(VanPrec, day.5, VanPrecfd, titles=CanadianWeather$place[26])

#  set up the functional parameter object for (positive smoothing

dayfdPar <- fdPar(smallbasis, harmaccelLfd365, lambda)

#  smooth the data with a positive function

Wfd1 <- smooth.pos(day.5, VanPrec, dayfdPar)$Wfdobj

#  plot both the original smooth and positive smooth

VanPrecvec    <- eval.fd(day.5, VanPrecfd)
VanPrecposvec <- eval.posfd(day.5, Wfd1)

par(ask=FALSE)
plot(day.5,  VanPrec, type="p")
lines(day.5, VanPrecposvec, lwd=2, col=4)
lines(day.5, VanPrecvec, lwd=2, col=3)
legend(100, 8, c("Positive smooth", "Unrestricted smooth"),
       lty=c(1,1), col=c(4,3))

#  plot the residuals

VanPrecres <- VanPrec - VanPrecposvec
plot(day.5, VanPrecres^2, type="p")
title("Squared residuals from positive fit")

#  compute a postive smooth of the squared residuals

lambda <- 1e3
dayfdPar <- fdPar(smallbasis, harmaccelLfd365, lambda)

Wfd <- smooth.pos(day.5, VanPrecres^2, dayfdPar)$Wfdobj

#  plot the square root of this smooth along with the residuals

VanPrecvarhat <- eval.posfd(day.5, Wfd)
VanPrecstdhat <- sqrt(VanPrecvarhat)

plot(day.5, VanPrecres^2, type="p")
lines(day.5, VanPrecvarhat, lwd=2)

#  set up a weight function for (revised smoothing

wtvec <- 1/VanPrecvarhat

lambda   <- 1e3
dayfdPar <- fdPar(smallbasis, harmaccelLfd365, lambda)

Wfd2 <- smooth.pos(day.5, VanPrec, dayfdPar, wtvec)$Wfdobj

#  plot the two smooths, one with weighting, one without

VanPrecposvec2 <- eval.posfd(day.5, Wfd2)

plot(day.5,  VanPrec, type="p")
lines(day.5, VanPrecposvec2, lwd=2, col=4)
lines(day.5, VanPrecposvec, lwd=2, col=3)
legend(100, 8, c("Weighted", "Unweighted"),
       lty=c(1,1), col=c(4,3))
##################################

     daybasis65 <- create.fourier.basis(c(0, 365), 65)

     daytempfd <- smooth.basis(day.5,
         CanadianWeather$dailyAv[,,"Temperature.C"],
         daybasis65, fdnames=list("Day", "Station", "Deg C"))$fd

     with(CanadianWeather, plotfit.fd(y=dailyAv[,,"Temperature.C"],
            argvals=day.5, daytempfd, index=1, titles=place) )
     with(CanadianWeather, plotfit.fd(y=dailyAv[,,"Temperature.C"],
            argvals=day.5, daytempfd, index=c(2, 35), titles=place) )
     with(CanadianWeather, plotfit.fd(y=dailyAv[,,"Temperature.C"],
            argvals=day.5, daytempfd, index=c(2, 35), titles=place,
               col=c("black", "blue") ) )
###################################




smallbasis  <- create.fourier.basis(c(0, 365), 65)

VanPrec <- with(CanadianWeather,
                dailyAv[, place='Vancouver', 'Precipitation.mm'])
harmaccelLfd365 <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))
dayfdPar <- fdPar(smallbasis, harmaccelLfd365, 1000)
Wfd1 <- smooth.pos(day.5, VanPrec, dayfdPar)$Wfdobj
VanPrecposvec <- eval.posfd(day.5, Wfd1)
VanPrecres <- VanPrec - VanPrecposvec

Wfd <- smooth.pos(day.5, VanPrecres^2, dayfdPar)$Wfdobj
VanPrecvarhat <- eval.posfd(day.5, Wfd)
wtvec <- wtcheck(length(VanPrecvarhat), 1/VanPrecvarhat)
Wfd2 <- smooth.pos(day.5, VanPrec, dayfdPar, wtvec)$Wfdobj
