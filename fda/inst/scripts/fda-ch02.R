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
###     of thier analyses but not always the data.
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
###
### ch. 2.  Tools for exploring functional data
###
###
library(fda)
##
## Section 2.1.  Introduction 
##
# No figures in this section.  

##
## Section 2.2.  Some notation  
##
# No figures in this section.  

##
## Section 2.3.  Summary statistics for functional data 
##
# pp. 22-23, Figure 2.1.
# Mean and standard deviation of pinch force
# after alignment / registration.

matplot(pinchtime, pinch, type="l")

# experiment with lambda to get a plot that looks
# reasonably smooth but not too smooth
pinch.fd <- smooth.basisPar(pinchtime, pinch)
pinch.fd <- smooth.basisPar(pinchtime, pinch, lambda=1e-6)

plot(pinch.fd$fd)
# The default lambda = 1/diff(pinchtime) produces a straight line
# lambda = 1e-8 undersmooths
# lambda = 1e-6 comes close to matching Figure 2.1
#   ... but with discrepancies
#   suggesting that the 'pinch' data used for Figure 2.1
#   may have differed by more than just an offset of 2
#   from the 'pinch' data in the 'fda' package 

plot(mean(pinch.fd$fd), xlab="Time (sec.)", ylab="Force mean (N)")
abline(v=.1, lty="dotted")

plot(sd(pinch.fd$fd), ylim=c(0, 1), xlab="Time (sec.)",
     ylab="Force Std. Dev. (N)")
abline(v=.1, lty="dotted")

# Reconsider after ch. 7?  

# pp. 23, Figure 2.2.  Perspective and contour plots of pinch force correlations

# complete after ch. 7 






# pp. 24-25, Figure 2.3.
# Contour plots of correlation and cross correlations of gait data

#tempPrecCor4.6 <- cor.fd(weeks, tempSmooth4$fd, fdobj2=precSmooth6$fd)

# In \demo\gaitDemo.R, lambda=1e-11 was determined as optimal smoothing 
# to minimize the generalized cross validation (GCV) criterion

# 'gaittime' is retained as the first of the gait dimnames 
(gaittime <- as.numeric(dimnames(gait)[[1]]))

# Establish a full basis;  
# 11 or 15 basis function might do as well or better,
# but I haven't tried that.  
gaitbasis21 <- create.fourier.basis(nbasis=21)

# A finite Fourier series must satisfy a harmonic
# differential equation, as described in the book.
# Give that a name.  
harmaccelLfd01 <- vec2Lfd(c(0, 0, (2*pi)^2, 0))

# lambda = 1e-11 is the roughness penalty
# found to minimize the generalized cross validation criterion;
# see '~R\library\demo\gaitDemo.R'.  
gaitfd <- smooth.basisPar(gaittime, gait,
       gaitbasis21, Lfdobj=harmaccelLfd01, lambda=1e-11)$fd

gaitfd$fdnames
# Change some of the names 
names(gaitfd$fdnames) = c("Normalized time", "Child", "Angle")
gaitfd$fdnames[[3]] = c("Hip", "Knee")

# Establish a finer grid so the plot appears smoother 
(gaitTimeSm <- seq(0, 1, length=41))
gaitCor <- cor.fd(gaitTimeSm, gaitfd)

gaitCorLbls <- (-5:5)/5
contour(gaitTimeSm, gaitTimeSm, gaitCor[,,1,1], cex=1.2,
        levels=gaitCorLbls, xlab="Time in Knee Angle Cycle",
        ylab="Time in Knee Angle Cycle")
x01. <- seq(0, 1, .05)
text(x01., x01., 1)
title("Knee - Knee")
# Smoothed more than Figure 2.3,
# which seems too rough 

contour(gaitTimeSm, gaitTimeSm, gaitCor[,,1,2], cex=1.2,
        levels=gaitCorLbls, ylab="Time in Knee Angle Cycle",
        xlab="Time in Hip Angle Cycle")
title("Hip - Knee")

contour(gaitTimeSm, gaitTimeSm, gaitCor[,,1,2], cex=1.2,
        levels=gaitCorLbls, xlab="Time in Knee Angle Cycle",
        ylab="Time in Hip Angle Cycle")
title("Knee - Hip")

contour(gaitTimeSm, gaitTimeSm, gaitCor[,,1,3], cex=1.2,
        levels=gaitCorLbls, xlab="Time in Hip Angle Cycle",
        ylab="Time in Hip Angle Cycle")
text(x01., x01., 1)
title("Hip - Hip")

##
## pp. 25-26, Figure 2.4.  Correlations for Canadian temperature & precipitation
##

# from demo\CanadianWeatherDemo.R
#
#  ----------------------------------------------------------------------
#                Descriptive Statistics Functions
#  ----------------------------------------------------------------------

#  -------------  set up fourier basis  ---------------------------
#  It was decided that 65 basis functions captured enough of 
#  the detail in the temperature data: about one basis function
#  per week for present purposes.

daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65)
#  -----------  set up the harmonic acceleration operator  ----------
harmaccelLfd365 <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))

# For other purposes, a saturated basis (365 basis functions) is used,
# where smoothing is defined by the GCV (Generalized Cross Validation) 
# criterion.

#  ---------  create fd objects for temp. and prec. ---------------

qqnorm(CanadianWeather$dailyAv[,,"Temperature.C"], datax=TRUE)
# Short tailed distribution
# because it's a strong, deterministic annual cycle
#   + relatively modest normal error

daytempfd <- with(CanadianWeather, data2fd(dailyAv[,,"Temperature.C"],
       day.5, daybasis65, argnames=list("Day", "Station", "Deg C")) )

qqnorm(CanadianWeather$dailyAv[,,"Precipitation.mm"], datax=TRUE)
# Apparent lognormal distribution

qqnorm(CanadianWeather$dailyAv[,,"log10precip"], datax=TRUE)
# Looks more like the 'tempav' distribution.
# Too discrete in the lower tail,
# because the numbers were rounded off too grossly
# before computing the logarithms.
sum(CanadianWeather$dailyAv[,,"Precipitation.mm"]==0)
# 27 of 365 numbers are 0 
quantile(CanadianWeather$dailyAv[,,"log10precip"])

# Use log10precip rather than Precipitation.mm directly.  

# Smooth temparature
# CanadianWeatherDemo.R compares alternative smoothing levels
# with a full basis rather than daybasis65.
# The generalized cross validation criterion is minimized for
# lambda = 0.01.  However, that ignores serial correlation and
# seems visually to undersmooth;  lambda = 1e6 looks more realistic
# and will be used here, even though the reduction in the number
# of basis vectors from 365 to 65 provides some smoothing without
# a roughness penalty.    
TempSm6 <- smooth.basisPar(argvals=day.5, y=CanadianWeather$dailyAv[,,1],
     fdobj=daybasis65, Lfdobj=harmaccelLfd365, lambda=1e6)$fd
# The minimal GCV in CanadianWeatherDemo.R for log10precip
# was also lambda=1e6 and will be used here as it gives a reasonable
# physical image in the plots.  
PrecipSm6 <- smooth.basisPar(argvals=day.5, y=CanadianWeather$dailyAv[,,3],
     fdobj=daybasis65, Lfdobj=harmaccelLfd365, lambda=1e6)$fd

lvl.1 <- (-10:10)/10
contour(weeks, weeks, cor.fd(weeks, TempSm6), levels=lvl.1,
        xlab="Average daily Temperature (C)", 
        ylab="Average daily Temperature (C)", 
        axes=FALSE, main="Temperature correlations")
axisIntervals(1)
axisIntervals(2)

contour(weeks, weeks, cor.fd(weeks, PrecipSm6), levels=lvl.1,
        xlab="Average daily Temperature (C)", 
        ylab="Average daily Temperature (C)", 
        axes=FALSE, main="Correlations in log10(Precipitation(mm))")
axisIntervals(1)
axisIntervals(2)

contour(weeks, weeks, cor.fd(weeks, TempSm6, weeks, PrecipSm6),
        levels=lvl.1, xlab="Average daily Temperature (C)", 
        ylab="Average daily Temperature (C)", axes=FALSE,
        main="Correlations in log10(Precipitation(mm))")
axisIntervals(1)
axisIntervals(2)
# Places where January is warmer are quite likely to have more rain (r=.8).
# Places where June is warmer are only moderately to have more rain (r=.4).


# This is roughly similar to Figure 2.4.
# To get the exact image,
# we would need the right combination of the following:  
# (1) The exact same basis set,
#     which probably was create.fourier.basis,
#     but might have been something else.
# (2) Exactly the same order 'nbasis'.
# (3) Exactly the same lambda and Lfdobj for smoothing.
# (4) Evalated at exactly the same points

# I don't know how to find that exact combination.
# However, I believe this is essentially equivalent
# to how that figure was created.  

##
## Section 2.4.  The anatomy of a function  
##

# No figures in this section

##
## Section 2.5.  Phase-plane plots of periodic effects 
##

# Figure 2.5.  Nondurable goods production in the US

plot(nondurables, log="y")
str(nondurables)

plot(log10(nondurables))
(log.nondur.lm <- lm(log(nondurables)~time(nondurables)))
(log10.nondur.lm <- lm(log10(nondurables)~time(nondurables)))
abline(nondur.lm, lty="dashed")

length(nondurables)

diff(range(time(nondurables)))
# 81

(nondurGrowth <- (nondurables[973]/nondurables[1]))
# 15.221 = (1+growthRate)^(81 years)
#log(1+growthRate) = log(15.221)/81 = 0.0336 
#    1+growthRate = nondurGrowth^(1/81) = 1.0342
#1.016^81 = 3.62

#1.0342^81 = 15.24

log(nondurables[973]/nondurables[1])/81

# Figure 2.6.  log nondurable goods 1964-1967 

nonDur1964.67 <- window(nondurables, 1964, 1967)

plot(log10(nondur1964.1967), type="p", axes=FALSE, xlab="Year",
     ylab=expression(paste(log[10], " nondurable goods index")) )
axis(2)
axis(1, 1964:1967)
axis(1, seq(1964, 1967, by=0.5), labels=FALSE)

# For theory on smoothing, see FDA, ch. 5
# For more on this example, see AFDA, ch. 3
# or \demo\goodindexDemo.R 

durtimefine <- seq(1964, 1967, length=181)

#fit = eval.fd(durtimefine, lognondursmth);
logNondurSm1964.67 = eval.fd(durtimefine, logNondurSm$fd);
lines(durtimefine, logNondurSm1964.67)
abline(v=1965:1966, lty=2)

# Figure 2.7.  Phase-plane plot for a simple harmonic function
sin. <- expression(sin(2*pi*x))
D.sin <- D(sin., "x")
D2.sin <- D(D.sin, "x")

with(data.frame(x=seq(0, 1, length=46)),
     plot(eval(D.sin), eval(D2.sin), type="l",
          xlim=c(-10, 10), ylim=c(-50, 50), 
          xlab="Velocity", ylab="Acceleration") )
pi.2 <- (2*pi)
#lines(x=c(-pi.2,pi.2), y=c(0,0), lty=3)
abline(h=0, lty="longdash")
pi.2.2 <- pi.2^2
lines(x=c(0,0), y=c(-pi.2.2, pi.2.2), lty="longdash")

text(c(0,0), c(-47, 47), rep("no kinetic, max potential", 2))
text(c(-8.5,8.5), c(0,0), rep("max kinetic\nno potential", 2))

# Figure 2.8.  Phase-plane plot of nondurable goods index 1964

goodsbasis <- create.bspline.basis(rangeval=c(1919,2000),
                                   nbasis=979, norder=8)
LfdobjNonDur <- int2Lfd(4) 

library(zoo)
logNondurSm <- smooth.basisPar(argvals=index(nondurables),
                y=log10(coredata(nondurables)), fdobj=goodsbasis,
                Lfdobj=LfdobjNonDur, lambda=1e-11)
phaseplanePlot(1964, logNondurSm$fd)

