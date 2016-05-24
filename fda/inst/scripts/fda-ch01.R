###
###
### Ramsey & Silverman (2006) Functional Data Analysis, 2nd ed. (Springer)
###
### DISCLAIMER:
###
### These script files are incomplete and where complete
### do not always reproduce exactly the images in the books.
### Reason for this include the following:
###
### (1) The authors were given permission to publish the results
###     of their analyses but not always the data.
### (2) These script files were produced starting in 2007
###     not by the authors to the books but by volunteer(s)
###     (initially Spencer Graves).  The published analyses
###     were conducted over more than a decade, and the details
###     of exactly how the published analyses were performed were
###     not in some cases still available.
### (3) The software and capabilities of computers available
###     have improved over time.  In some cases, this produces
###     a better fit today than at the times the published analyses
###     were prepared.
###
### It is hoped that these script files may still be useful,
### even with their gaps and even when they do not reproduce exactly
### the published analyses.
###
### Spencer Graves
###
### ch. 1.  Introduction
###
###
library(fda)
##
## Section 1.1.  What are functional data?
##
# pp. 1-2, Figure 1.1.  Heights of 10 girls measured at 31 ages
#str(growth)
#sapply(growth, class)
with(growth, matplot(age, hgtf[, 1:10], type="b"))

# Compare on different scales:
op <- par(mfrow=c(2,2))

# connecting the dots
with(growth, matplot(age, hgtf[, 1:10], type="b",
                     main='linear'))
#with(growth, matplot(age, hgtf, type="b"))

with(growth, matplot(age, hgtf[, 1:10], type="b", log='x',
                     main="x log scale"))
with(growth, matplot(age, hgtf[, 1:10], type="b", log='y',
                     main="y log scale"))
with(growth, matplot(age, hgtf[, 1:10], type="b", log='xy',
                     main="log-log"))
par(op)
# Best scale (most nearly linear) of these 4:  log-log

# Alternative:  with smoothing
# Use bspline basis with a knot at each age
girlGrowthSm <- with(growth, smooth.basisPar(age, hgtf))

sapply(girlGrowthSm, class)
plot(girlGrowthSm$fd[1:10])
with(growth, matpoints(age, hgtf[, 1:10]))

# alternative plot

with(growth, plotfit.fd(hgtf[, 1:10], age, girlGrowthSm$fd[1:10]))

# Matches Figure 1.1 quite well.


# Figure 1.2.  Estimated accelerations for the heights of 10 girls

plot(deriv(girlGrowthSm$fd[1:10], 2))
# NOT smooth enough

# Experiment with the level of smoothing via lambda:
girlGrowthSm.25 <- with(growth, smooth.basisPar(age, hgtf, lambda=.25))
plot(deriv(girlGrowthSm.25$fd[1:10], 2), ylim=c(-4, 2))
# General match to Figure 1.2 but off in some details.
# To get a better match, we may need to use a different smoothing
# operator in place of the default 'Lfdobj=int2Lfd(2)'.
lines(mean(deriv(girlGrowthSm.25$fd[1:10], 2)), lty="dashed", lwd=3)

# pp. 3-4, Figure 1.3.  US nondurable goods manufacturing index
Nondurables <- ts(nondurables, start=c(1919, 10), frequency=12)
plot(Nondurables, ylab="Nondurable Goods Index")

# Or on the log scale
plot(Nondurables, ylab="Nondurable Goods Index", log='y')

# Figure 1.4.  Tray level & vapor flow during a refinery experiment
op <- par(mfrow=c(2, 1), mar=c(2, 4, 4, 5)+0.1, cex.lab=1.5,
          las=1)
with(refinery, plot(Time, Tray47, pch=".", cex=2,
                    ylab="Tray 47 level"))
axis(3, labels=FALSE)
par(mar=c(5, 4, 1, 5)+0.1)
with(refinery, plot(Time, Reflux, pch=".", cex=3, xlab="Time",
                    ylab="Refulx flow"))
axis(3, labels=FALSE)
par(op)

##
## 1.2.  Functional models for nonfunctional data
##
# p. 5-6, Figure 1.5.  Item response vs. latent math ability

# Data not available.

##
## 1.3.  Some functional data analyses
##
# p. 5-6, Figure 1.6.  Mean monthly temperatures
#              for Canadian weather stations

# Approximate each curve of 365 points with a fit to
# a finite Fourier series of length 65
# = mean + 32 (sine-cosine) pairs
daybasis65 <- create.fourier.basis(rangeval=c(0, 365), nbasis=65)

# Smooth each curve additionally using a roughness penalty
# based on "Harmonic Acceleration";  see later in FDA:
harmaccelLfd365 <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))

# Use a roughness penalty of 1e6;
# see \demo\CanadianWeatherDemo.R for a discussion of alternative
# smoothing penalties ...
# A penalty of 1e5 presents an image rougher than that in Figure 1.6.

CanadianWeather$place
Stations <- c(Pr.Rupert=29, Montreal=12, Edmonton=23, Resolute=35)
#attach(CanadianWeather)
TempSmooth6.4places <- smooth.basisPar(argvals=day.5,
     y=CanadianWeather$dailyAv[,Stations, "Temperature.C"],
     fdobj=daybasis65, Lfdobj=harmaccelLfd365, lambda=1e6)$fd

op <- par(mar=c(5, 4, 4, 5)+.1)
plot(TempSmooth6.4places, axes=FALSE,
     xlab="Month", ylab="Temperature (C)", bty="n")
matpoints(monthMid, CanadianWeather$monthlyTemp[, Stations])
# lambda=1e5 presents a rougher image than Figure 1.6.

axisIntervals(1)
axis(2, las=1)
#axis(1, Month, substring(names(Month), 1,1))
text(rep(399, length(Stations)),
     CanadianWeather$dailyAv[365,Stations, "Temperature.C"],
     CanadianWeather$place[Stations], xpd=NA)

par(op)

# p. 7, Figure 1.7.  The results of applying
# the differential operator L = (pi/6)^2*D+D^3
# to the estimated temperature functions in Figure 1.6

plot(TempSmooth6.4places, Lfdobj=harmaccelLfd365,
     axes=FALSE, xlab="Month", ylab="L-Temperature", bty="n")

axisIntervals(1)
axis(2, las=1)
#axis(1, Month, substring(names(Month), 1,1))

#The results of applying the harmonic differential operator
# L = (pi/6)^2*D+D^3 to the estimated temperature functions in Figure 1.6

# For Figure 1.7, ylim is approximately c(-22, 22).
# This script produces a plot with ylim approximately c(-5e-4, 5e-4).
# Either there is an error in Figure 1.7
# or Figure 1.6 was not produced using the harmonic differential operator.


# p. 8, Figure 1.8.  Hip and knee angles through a gait cycle

str(gait)
op <- par(mfrow=c(2, 1), mar=c(2, 4, 4, 5)+0.1, cex.lab=1.5,
          las=1)
# Hip Angle
plot(c(-.2, 1.2), range(gait[,,"Hip Angle"]), type = "n",
     ylab="Hip angle (degrees)")

ind1 <- c(20, 1:20)
matlines((0:20)/20, gait[ind1, , "Hip Angle"], lty=1)
ind0 <- (-(4:0))
matlines(ind0/20, gait[20+ind0, , "Hip Angle"], lty=3)
ind2 <- c(20, 1:4)
matlines((20:24)/20, gait[ind2, , "Hip Angle"], lty=3)
abline(v=0:1, lty=3)
# Knee Angle
par(mar=c(5, 4, 1, 5)+0.1)
plot(c(-.2, 1.2), range(gait[,,"Knee Angle"]), type = "n",
     ylab="Knee angle (degrees)",
     xlab="Time (proportion of gait cycle)")
matlines((0:20)/20, gait[ind1, , "Knee Angle"], lty=1)
ind0 <- (-(4:0))
matlines(ind0/20, gait[20+ind0, , "Knee Angle"], lty=3)
ind2 <- c(20, 1:4)
matlines((20:24)/20, gait[ind2, , "Knee Angle"], lty=3)
abline(v=0:1, lty=3)

par(op)

##
## 1.4.  The goals of a functional analysis
##
# p. 9-10.  Figure 1.9.  Gait cycle in the sagittal plane

plot(range(gait[,,"Hip Angle"]), range(gait[,,"Knee Angle"]),
     xlab="Hip angle (degrees)", ylab="Knee angle (degrees)",
     type="n")
subj <- 1
lines(gait[,subj,], type="b")
ltrs <- 4*(1:5)
text(gait[ltrs,subj,], LETTERS[c(2:5, 1)])

(gaitMean <- apply(gait, c(1, 3), mean) )
lines(gaitMean, lty=3, type="b")
text(gaitMean[ltrs,], LETTERS[c(2:5, 1)])

##
## 1.5.  The  first steps in a functional data analysis
##
# p. 11-12, Figure 1.10.  Average daily rainfall at Prince Rupert

str(CanadianWeather)
CanadianWeather$place
Pr.Rupert <- 29
#attach(CanadianWeather)
daybasis65 <- create.fourier.basis(rangeval=c(0, 365), nbasis=65)
harmaccelLfd365 <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))

precipSm7.Pr.Rupert <- smooth.basisPar(argvals=day.5,
     y=CanadianWeather$dailyAv[,Pr.Rupert, "Precipitation.mm"],
     fdobj=daybasis65, Lfdobj=harmaccelLfd365, lambda=1e7)$fd

plotfit.fd(CanadianWeather$dailyAv[,Pr.Rupert, "Precipitation.mm"],
           day.5, precipSm7.Pr.Rupert, axes=FALSE, bty="n",
           ylab="Precipitation.mm", xlab="day of the year",
           main="Average daily rainfall at Prince Rupert")

axisIntervals(1)
axis(2, las=1)
#axis(1, Month, substring(names(Month), 1,1))


# p. 12-13, Figure 1.11.  Force exerted by thumb and forefinger
# Data not available.

# 'pinch' = (x-2), adjusted to a common start time
# from what is plotted in Figure 1.11.

matplot(pinchtime, pinch, xlab="Seconds", ylab="Force (N)", type="l")


# p. 14.  Figure 1.12.  Annual variation in mean temperature at Montreal

#op <- par(mfrow=c(1, 2), cex=2)
op <- par(cex=2)
#str(CanadianWeather)
tempMontreal <- CanadianWeather$dailyAv[, "Montreal", "Temperature.C"]

tempSm6.Montreal <- smooth.basisPar(argvals=day.5,
     y=CanadianWeather$dailyAv[,"Montreal", "Temperature.C"],
     fdobj=daybasis65, Lfdobj=harmaccelLfd365, lambda=1e6)$fd

plot(tempSm6.Montreal, axes=FALSE,
     xlab="Month", ylab="Mean Temperature")
axisIntervals(1)
axis(2)
text(monthMid, CanadianWeather$monthlyTemp[, "Montreal"],     monthLetters)
title(paste(c("Montreal average daily temp\ndeviation from average",
              " (C)"), collapse=""))

phaseplanePlot(fdobj=tempSm6.Montreal, Lfdobj1=0,
               xlab="Temperature (C)")
title(paste(c("Montreal average daily temp\ndeviation from average",
              " (C)"), collapse=""))
par(op)

# NOTE:  Figure 1.12 has the wrong sign on one of the axes:
# The slope should be negative, since
# D^2(sin) = -sin & D^2(cos) = -cos ...


# p. 14-15.  Figure 1.13.
# Phase plan plot of first 2 derivatives of US nondurable goods index, 1964
goodsbasis <- create.bspline.basis(rangeval=c(1919,2000),
                                   nbasis=979, norder=8)
LfdobjNonDur <- int2Lfd(4)

library(zoo)
logNondurSm <- smooth.basisPar(argvals=index(nondurables),
                y=log10(coredata(nondurables)), fdobj=goodsbasis,
                Lfdobj=LfdobjNonDur, lambda=1e-11)
phaseplanePlot(1964, logNondurSm$fd)

