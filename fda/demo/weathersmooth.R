
#  This file is intended to be used after the weather data have been set up
#  by the commands in weathersetup.R

#  The weather data are smoothed using the harmonic acceleration penaltym
#  and the smoothed data are saved as a named list weatherfd.

#  Many other interesting ways of smoothing the data and plotting results
#  can be found in the file canadian-weather.R, set up in 2008 by 
#  Spencer Graves.

#  Last modified 17 November 2008

#  load the data

#load("weatherdata")

#  ------------------------  set up the data  -----------------------

tempav  <- daily$tempav
precav  <- daily$precav
station <- daily$place

#  set up the times of observation at noon

daytime   <- (1:365)-0.5

daytime <- weatherdata$daytime
station <- weatherdata$station

#  JJindex re-orders days so that they run from July 1 to June 30.
#  This can be useful for understanding weather data because variation 
#  from curve to curve is greater and more interesting in mid-winter.  
#  July 1 is the 182nd day.
#  If this reordering is desired, use commands
#  daytime = daytime[JJindex]
#  tempav = tempav[JJindex,]
#  precav = precav[JJindex,]

JJindex = c(182:365, 1:181)

#  -------------  set up fourier basis  ---------------------------

#  Here it was decided that 65 basis functions captured enough of
#  the detail in the temperature data: about one basis function
#  per week.  However, see below for smoothing with a saturated
#  basis (365 basis functions) where smoothing is defined by the
#  GCV criterion.

rangeval <- c(0,365)
period   <- 365
nbasis   <- 365
daybasis <- create.fourier.basis(rangeval, nbasis, period)

#  -----------  set up the harmonic acceleration operator  ----------

harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), rangeval)

#  ---------------------------------------------------------------------
#                        Smooth temperature  
#  ---------------------------------------------------------------------

tempav <- weatherdata$tempav

#                 Choose level of smoothing using
#          the generalized cross-validation criterion
#              with smoothing function smooth.basis.

#  set up range of smoothing parameters in log_10 units

loglam <- -3:3
nlam   <- length(loglam)

dfsave  <- rep(0,nlam)
gcvsave <- rep(0,nlam)

#  loop through smoothing parameters

for (ilam in 1:nlam) {
	lambda     <- 10^loglam[ilam]
	cat(paste("lambda =",lambda,"\n"))
	fdParobj   <- fdPar(daybasis, harmaccelLfd, lambda)
	smoothlist <- smooth.basis(daytime, tempav, fdParobj)
	fdobj      <- smoothlist[[1]]
	df         <- smoothlist[[2]]
	gcv        <- smoothlist[[3]]
	dfsave[ilam]  <- df
	gcvsave[ilam] <- sum(gcv)
}

#  display and plot degrees of freedom and GCV criterion

cbind(loglam, dfsave, gcvsave)

par(mfrow=c(1,2), pty="m")
plot(loglam, gcvsave, type="b", cex=1,
     xlab="Log_10 lambda", ylab="GCV Criterion",
     main="Temperature Smoothing")

plot(loglam, dfsave, type="b",  cex=1,
     xlab="Log_10 lambda", ylab="Degrees of freedom",
     main="Temperature Smoothing")

#  Do final smooth with lambda = 100 because the minimum GCV value
#  is very rough, and also equivalent to 255 degrees of freedom,
#  whereas the smooth with this smoothing parameter value is able 
#  to nicely track high frequency events such as the January thaws
#  clearly visible in the data.  Moreover, the lambda = 100 smooth
#  is equivalent to about 56 degrees of freedom, which is enough to
#  track events down to a resolution about one week.

lambda   <- 100  #  minimum GCV estimate, corresponding to 255 df
fdParobj <- fdPar(daybasis, harmaccelLfd, lambda)

smoothlist <- smooth.basis(daytime, tempav, fdParobj)

daytempfd  <- smoothlist$fd
df         <- smoothlist$df
SSE        <- smoothlist$SSE
tempy2cMap <- smoothlist$y2cMap

print(paste("Degrees of freedom =", round(df,0)))

daytempfd$fdnames <- list(Day = NULL, Station = station, DegC = NULL)

#  estimate standard error of fit

stderr <- sqrt(SSE/(35*(365-df)))  

print(paste("Standard error =", round(stderr,2), "deg. C"))

#  plot data and fit

par(mfrow=c(1,1), pty="m", ask=T)
plotfit.fd(tempav, daytime, daytempfd, titles=station)

#  Compute the standard error of estimate for the fitted curve

Phimat <- eval.basis(daytime, daybasis)

yhatMap <- Phimat %*% tempy2cMap

SigmaErrMat <- crossprod(t(yhatMap))*stderr^2

StdErrVec <- sqrt(diag(SigmaErrMat))

#  Set up the data for the winter days, from Dec. 15 to Feb. 28

WinterIndex <- c(340:365,1:66)
WinterTime  <- (1:length(WinterIndex))-0.5

#  Plot the winter data, the fit, and 
#  the pointwise 95% confidence intervals for the fit

tempavhat <- eval.fd(daytime, daytempfd)
par(ask=T)
for (i in 1:35) {
  plot(WinterTime, tempavhat[WinterIndex,i], type="l", cex = 1.2,
       xlim=c(0,length(WinterIndex)), ylim=c(-35,5),
       xlab="Day from Dec. 15", ylab="Deg. C", main=station[i])
  lines(WinterTime, tempavhat[WinterIndex,i]+2*StdErrVec[WinterIndex], lty=3)
  lines(WinterTime, tempavhat[WinterIndex,i]-2*StdErrVec[WinterIndex], lty=3)
  points(WinterTime, tempav[WinterIndex,i])
}
      
#  ---------------------------------------------------------------------
#                        Smooth precipitation  
#  ---------------------------------------------------------------------

#                 Choose level of smoothing using
#          the generalized cross-validation criterion
#              with smoothing function smooth.basis.

#  set up range of smoothing parameters in log.10 units

loglam <- 4:9
nlam <- length(loglam)

dfsave  <- rep(0,nlam)
gcvsave <- rep(0,nlam)

#  loop through smoothing parameters

for (ilam in 1:nlam) {
	lambda <- 10^loglam[ilam]
	cat(paste("lambda =",lambda,"\n"))
	fdParobj <- fdPar(daybasis, harmaccelLfd, lambda)
	smoothlist <- smooth.basis(daytime, precav, fdParobj)
	fdobj  <- smoothlist[[1]]
	df     <- smoothlist[[2]]
	gcv    <- smoothlist[[3]]
	dfsave[ilam]  <- df
	gcvsave[ilam] <- sum(gcv)
}

#  display and plot degrees of freedom and GCV criterion

cbind(loglam, dfsave, gcvsave)

par(mfrow=c(1,2), pty="m")
plot(loglam, gcvsave, type="b", cex=1,
     xlab="Log_10 lambda", ylab="GCV Criterion",
     main="Precipitation Smoothing")

plot(loglam, dfsave, type="b", cex=1,
     xlab="Log_10 lambda", ylab="Degrees of freedom",
     main="Precipitation Smoothing")

#  Do final smooth with lambda = 1e5 rather than minimum GCV value
#  of 1e7 because the latter seems to over-smooth several sets of data

lambda   <- 1e5  
fdParobj <- fdPar(daybasis, harmaccelLfd, lambda)

smoothlist <- smooth.basis(daytime, precav, fdParobj)

dayprecfd  <- smoothlist$fd
df         <- smoothlist$df
SSE        <- smoothlist$SSE
precy2cMap <- smoothlist$y2cMap

print(paste("Degrees of freedom =", round(df,0)))

dayprecfd$fdnames <- list(Day = NULL, Station = station, Millimetres = NULL)

#  estimate standard error of fit

stderr <- sqrt(SSE/(35*(365-df))) 

print(paste("Standard error =", round(stderr,2), "mm"))

#  plot data and fit

par(mfrow=c(1,1), pty="m", ask=T)
plotfit.fd(precav, daytime, dayprecfd, titles=station)

#  ---------------------------------------------------------------------
#    Plot precipitation against temperature for each station in turn
#  ---------------------------------------------------------------------

daytime <- weatherdata$daytime
station <- weatherdata$station
tempmat <- eval.fd(daytime, daytempfd)
precmat <- eval.fd(daytime, dayprecfd)

#  define 1-character names for months

monthletter <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
monthtime   <- c(1,32,50,81,111,142,172,203,234,264,295,315)

par(ask=T)
for (i in 1:35) {
  plot(tempmat[,i], precmat[,i], type="l", cex=1.2, 
       xlim=c(-35,25), ylim=c(0,13), 
       xlab="Temperature", ylab="Precipitation",
       main=station[i])
  text(tempmat[monthtime,i], precmat[monthtime,i], monthletter, cex=1.2)
}

#  -------------------------------------------------------------------
#         save the results as two fdPar objects plus y2cMaps
#  -------------------------------------------------------------------

tempfdPar <- fdPar(daytempfd, harmaccelLfd, 1e2)
precfdPar <- fdPar(dayprecfd, harmaccelLfd, 1e5)

weatherfd <- list(tempfdPar  = tempfdPar,  precfdPar  = precfdPar,
                  tempy2cMap = tempy2cMap, precy2cMap = precy2cMap)

save(weatherfd, file="weatherfd")

#  -------------------------------------------------------------------
#     Smooth Vancouver's precipitation with a positive function.
#  -------------------------------------------------------------------

#  select Vancouver's precipitation

index <- (1:35)[station == "Vancouver  "] 
VanPrec  <- precav[,index]

#  We use a more compact basis here to avoid storage allocation
#  problems that arise when 365 basis functions are used.

dayrange   <- c(0,365)
daybasis65 <- create.fourier.basis(dayrange, 65)

#  smooth the data using 65 basis functions

lambda    <- 1e4
fdParobj  <- fdPar(daybasis65, harmaccelLfd, lambda)
VanPrecfd <- smooth.basis(daytime, VanPrec, fdParobj)$fd

#  Plot temperature curves and values

plotfit.fd(VanPrec, daytime, VanPrecfd, titles=station[index])

#  smooth the data with a positive function

Wfd1 <- smooth.pos(daytime, VanPrec, fdParobj)$Wfdobj

#  plot both the original smooth and positive smooth

VanPrecvec    <- eval.fd(daytime, VanPrecfd)
VanPrecposvec <- eval.posfd(daytime, Wfd1)

par(ask=F)
plot(daytime,  VanPrec, type="p",
     xlab="Day", ylab="Precipitation (mm)")
lines(daytime, VanPrecposvec, lwd=2, col=4)
lines(daytime, VanPrecvec, lwd=2, col=2)
legend(100, 8, c("Positive smooth", "Unrestricted smooth"), 
       lty=c(1,1), col=c(4,2))

#  plot the squared residuals

VanPrecres <- VanPrec - VanPrecposvec
plot(daytime, VanPrecres^2, type="p",
     xlab="Day", ylab="Squared Residuals")
title("Squared residuals from positive fit")

#  compute a postive smooth of the squared residuals

lambda <- 1e3
fdParobj <- fdPar(daybasis65, harmaccelLfd, lambda)

Wfd <- smooth.pos(daytime, VanPrecres^2, fdParobj)$Wfdobj

#  plot the square root of this smooth along with the residuals

VanPrecvarhat <- eval.posfd(daytime, Wfd)
VanPrecstdhat <- sqrt(VanPrecvarhat)

plot(daytime, VanPrecres^2, type="p",
     xlab="Day", ylab="Squared Residuals")
lines(daytime, VanPrecvarhat, lwd=2)

#  set up a weight function for (revised smoothing

wtvec <- as.vector(1/VanPrecvarhat)

lambda   <- 1e3
fdParobj <- fdPar(daybasis65, harmaccelLfd, lambda)

Wfd2 <- smooth.pos(daytime, VanPrec, fdParobj, wtvec)$Wfdobj

#  plot the two smooths, one with weighting, one without

VanPrecposvec2 <- eval.posfd(daytime, Wfd2)

plot(daytime,  VanPrec, type="p")
lines(daytime, VanPrecposvec2, lwd=2, col=4)
lines(daytime, VanPrecposvec, lwd=2, col=2)
legend(100, 8, c("Weighted", "Unweighted"), 
       lty=c(1,1), col=c(4,2))



