###
### Ramsay, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###

#  Remarks and disclaimers

#  These R commands are either those in this book, or designed to
#  otherwise illustrate how R can be used in the analysis of functional
#  data.
#  We do not claim to reproduce the results in the book exactly by these
#  commands for various reasons, including:
#    -- the analyses used to produce the book may not have been
#       entirely correct, possibly due to coding and accuracy issues
#       in the functions themselves
#    -- we may have changed our minds about how these analyses should be
#       done since, and we want to suggest better ways
#    -- the R language changes with each release of the base system, and
#       certainly the functional data analysis functions change as well
#    -- we might choose to offer new analyses from time to time by
#       augmenting those in the book
#    -- many illustrations in the book were produced using Matlab, which
#       inevitably can imply slightly different results and graphical
#       displays
#    -- we may have changed our minds about variable names.  For example,
#       we now prefer "yearRng" to "yearRng" for the weather data.
#    -- three of us wrote the book, and the person preparing these scripts
#       might not be the person who wrote the text
#  Moreover, we expect to augment and modify these command scripts from time
#  to time as we get new data illustrating new things, add functionality
#  to the package, or just for fun.

###
### ch. 4  How to Build Functional Data Objects
###

#  load the fda package

library(fda)

#  display the data files associated with the fda package

data(package='fda')

#  start the HTML help system if you are connected to the Internet, in
#  order to open the R-Project documentation index page in order to obtain
#  information about R or the fda package.

help.start()

##
## Section 4.1 Adding Coefficients to Bases to Define Functions
##

#  4.1.1 Coefficient Vectors, Matrices and Arrays

daybasis65 = create.fourier.basis(c(0,365), 65)
# dummy coefmat
coefmat = matrix(0, 65, 35, dimnames=list(
     daybasis65$names, CanadianWeather$place) )
tempfd. = fd(coefmat, daybasis65)

# 4.1.2 Labels for Functional Data Objects

fdnames      = list("Age (years)", "Child", "Height (cm)")

# or

fdnames      = vector('list', 3)
fdnames[[1]] = "Age (years)"
fdnames[[2]] = "Child"
fdnames[[3]] = "Height (cm)"

station      = vector('list', 35)
station[[ 1]]= "St. Johns"
#.
#.
#.
station[[35]] = "Resolute"

# Or:

station = as.list(CanadianWeather$place)

fdnames = list("Day", "Weather Station" = station,
    "Mean temperature (deg C)")

##
## 4.2 Methods for Functional Data Objects
##

#  Two order 2 splines over unit interval

unitRng = c(0,1)
bspl2 = create.bspline.basis(unitRng, norder=2)
plot(bspl2, lwd=2)

#  a pair of straight lines

tstFn1 = fd(c(-1, 2), bspl2)
tstFn2 = fd(c( 1, 3), bspl2)

#  sum of these straight lines

par(mfrow=c(3,1))
fdsumobj = tstFn1+tstFn2
plot(tstFn1,   lwd=2, xlab="", ylab="Line 1")
plot(tstFn2,   lwd=2, xlab="", ylab="Line 2")
plot(fdsumobj, lwd=2, xlab="", ylab="Line1 + Line 2")

#  difference between these lines

fddifobj = tstFn2-tstFn1
plot(tstFn1,   lwd=2, xlab="", ylab="Line 1")
plot(tstFn2,   lwd=2, xlab="", ylab="Line 2")
plot(fddifobj, lwd=2, xlab="", ylab="Line2 - Line 1")

fdprdobj = tstFn1 * tstFn2
plot(tstFn1,   lwd=2, xlab="", ylab="Line 1")
plot(tstFn2,   lwd=2, xlab="", ylab="Line 2")
plot(fdprdobj, lwd=2, xlab="", ylab="Line2 * Line 1")

#  square of a straight line

fdsqrobj = tstFn1^2
plot(tstFn1,   lwd=2, xlab="", ylab="Line 1")
plot(tstFn2,   lwd=2, xlab="", ylab="Line 2")
plot(fdsqrobj, lwd=2, xlab="", ylab="Line 1 ^2")

#  square root of a line with negative values:  illegal

a = 0.5
# fdrootobj = tstFn1^a
#Error in `^.fd`(tstFn1, a) :
#  There are negative values and the power is a positive fraction.

#  square root of a square:  this illustrates the hazards of
#  fractional powers when values are near zero.  The right answer is
#  two straight line segments with a discontinuity in the first
#  derivative.  It would be better to use order two splines and
#  put a knot at the point of discontinuity, but the power method
#  doesn't know how to do this.

fdrootobj = fdsqrobj^a
plot(tstFn1,    lwd=2, xlab="", ylab="Line 1")
plot(tstFn2,    lwd=2, xlab="", ylab="Line 2")
plot(fdrootobj, lwd=2, xlab="", ylab="sqrt(fdsqrobj)")

#  square root of a quadratic without values near zero:  no problem

fdrootobj = (fdsqrobj + 1)^a
par(mfrow=c(2,1))
plot(fdsqrobj + 1,    lwd=2, xlab="", ylab="fdsqrobj + 1")
plot(fdrootobj, lwd=2, xlab="", ylab="sqrt(fdsqrobj + 1)")

#  reciprocal of a function with zero values:  illegal operation

a    = (-1)
# fdinvobj = tstFn1^a
# Error in `^.fd`(tstFn1, a) :
#   There are zero or negative values and the power is negative.

#  reciprocal of a function with near zero values:  a foolish thing
#  to do and the power function fails miserably

fdinvobj = fdsqrobj^a
plot(fdsqrobj, lwd=2, xlab="", ylab="fdsqrobj")
plot(fdinvobj, lwd=2, xlab="", ylab="1/fdsqrobj")

#  reciprocal of a positive function with no values near zero

fdinvobj = (fdsqrobj+1)^a
plot(fdsqrobj + 1, lwd=2, xlab="", ylab="fdsqrobj + 1")
plot(fdinvobj,     lwd=2, xlab="", ylab="1/(fdsqrobj+1)")

#  near reciprocal of a positive function with no values near zero

a = -0.99
fdpowobj = (fdsqrobj+1)^a
plot(fdsqrobj + 1, lwd=2, xlab="", ylab="fdsqrobj + 1")
plot(fdpowobj,     lwd=2, xlab="", ylab="(fdsqrobj+1)^(-0.99)")

#
# compute mean temperature in two ways and plot the difference
#

Tempbasis = create.fourier.basis(yearRng, 65)
Tempfd = smooth.basis(day.5,
          CanadianWeather$dailyAv[,,'Temperature.C'], Tempbasis)$fd
meanTempfd = mean(Tempfd)
sumTempfd  = sum(Tempfd)

par(mfrow=c(1,1))
plot((meanTempfd-sumTempfd*(1/35)))

# round off error, as it should be.

#  plot the temperature for Resolute and add the Canadian mean

plot(Tempfd[35], lwd=1, ylim=c(-35,20))
lines(meanTempfd, lty=2)

#  evaluate the derivative of mean temperature and plot

DmeanTempVec = eval.fd(day.5, meanTempfd, 1)
plot(day.5, DmeanTempVec, type='l')

#  evaluate and plot the harmonic acceleration of mean temperature

harmaccelLfd = vec2Lfd(c(0,c(2*pi/365)^2, 0), c(0, 365))
LmeanTempVec = eval.fd(day.5, meanTempfd, harmaccelLfd)

par(mfrow=c(1,1))
plot(day.5, LmeanTempVec, type="l", cex=1.2,
     xlab="Day", ylab="Harmonic Acceleration")
abline(h=0)

#  plot Figure 4.1

dayOfYearShifted = c(182:365, 1:181)

tempmat   = daily$tempav[dayOfYearShifted, ]
tempbasis = create.fourier.basis(yearRng,65)

temp.fd = smooth.basis(day.5, tempmat, tempbasis)$fd

temp.fd$fdnames = list("Day (July 2 to June 30)",
                      "Weather Station",
                      "Mean temperature (deg. C)")

plot(temp.fd, lwd=2, xlab='Day (July 1 to June 30)',
     ylab='Mean temperature (deg. C)')

#
# Section 4.2.1 Illustration: Sinusoidal Coefficients
#

# Figure 4.2

basis13  = create.bspline.basis(c(0,10), 13)
tvec     = seq(0,1,len=13)
sinecoef = sin(2*pi*tvec)
sinefd   = fd(sinecoef, basis13, list("t","","f(t)"))
op       = par(cex=1.2)
plot(sinefd, lwd=2)
points(tvec*10, sinecoef, lwd=2)
par(op)

##
## Section 4.3 Smoothing using Regression Analysis
##

# Section 4.3.1 Plotting the January Thaw

# Figure 4.3

# This assumes the data are in "MtlDaily.txt"
# in the working directory getwd();
# first create it and put it there
cat(MontrealTemp, file='MtlDaily.txt')

MtlDaily = matrix(scan("MtlDaily.txt",0),34,365)
thawdata = t(MtlDaily[,16:47])

daytime  = ((16:47)+0.5)
plot(daytime, apply(thawdata,1,mean), "b", lwd=2,
     xlab="Day", ylab="Temperature (deg C)", cex=1.2)

# Figure 4.4

thawbasis    = create.bspline.basis(c(16,48),7)
thawbasismat = eval.basis(thawbasis, daytime)

thawcoef = solve(crossprod(thawbasismat),
    crossprod(thawbasismat,thawdata))
thawfd   = fd(thawcoef, thawbasis,
    list("Day", "Year", "Temperature (deg C)"))
plot(thawfd, lty=1, lwd=2, col=1)

# Figure 4.5

plotfit.fd(thawdata[,1], daytime, thawfd[1],
           lty=1, lwd=2, main='')

##
## Section 4.4 The Linear Differential Operator or Lfd Class
##

omega           = 2*pi/365
thawconst.basis = create.constant.basis(thawbasis$rangeval)

betalist       = vector("list", 3)
betalist[[1]]  = fd(0, thawconst.basis)
betalist[[2]]  = fd(omega^2, thawconst.basis)
betalist[[3]]  = fd(0, thawconst.basis)
harmaccelLfd.  = Lfd(3, betalist)

accelLfd = int2Lfd(2)

harmaccelLfd.thaw = vec2Lfd(c(0,omega^2,0), thawbasis$rangeval)
all.equal(harmaccelLfd.[-1], harmaccelLfd.thaw[-1])

class(accelLfd)
class(harmaccelLfd)

Ltempmat  = eval.fd(day.5, temp.fd, harmaccelLfd)

D2tempfd = deriv.fd(temp.fd, 2)
Ltempfd  = deriv.fd(temp.fd, harmaccelLfd)

##
## Section 4.5 Bivariate Functional Data Objects:
##             Functions of Two Arguments
##

Bspl2 = create.bspline.basis(nbasis=2, norder=1)
Bspl3 = create.bspline.basis(nbasis=3, norder=2)

corrmat  = array(1:6/6, dim=2:3)
bBspl2.3 = bifd(corrmat, Bspl2, Bspl3)

##
## Section4.6 The structure of the fd and Lfd Classes
##

help(fd)
help(Lfd)

##
## Section 4.7 Some Things to Try
##
# (exercises for the reader)
