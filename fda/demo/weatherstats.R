#  This file is intended to be used after the commands in files
#  weathersetup.R and weathersmooth.R have been executed.

#  Many other interesting ways of describing the data and plotting results
#  can be found in the file canadian-weather.R, set up in 2008 by 
#  Spencer Graves.

#  Compute and display various descriptive statistics for daily
#  temperature and precipitation

#  Last modified 17 November 2008

load("weatherfd")
load("weatherdata")

#  ---------  create fd objects for temp. and prec. ---------------

daytime   <- weatherdata$daytime
daytempfd <- weatherfd$tempfdPar$fd
dayprecfd <- weatherfd$precfdPar$fd

#  --  compute and plot mean and standard deviation of temperature -------

tempmeanfd <- mean(daytempfd)
tempstdvfd <- std.fd(daytempfd)

par(mfrow=c(2,1), pty="m", ask=F)
plot(tempmeanfd, main="Mean")
plot(tempstdvfd, main="Standard Deviation")

#  ---------  set up 53 day values roughly in weeks ----------------

weeks <- seq(0,365,length=53)

#  ---------------------------------------------------------------------
#    Plot by contour and perspective the variance-covariance and
#  correlation surfaces for both variables and for their cross-variation
#  ---------------------------------------------------------------------

#  --  plot the temperature variance-covariance bivariate function  ----

tempvarbifd <- var.fd(daytempfd)
tempvarmat  <- eval.bifd(weeks,weeks,tempvarbifd)

par(mfrow=c(1,1), pty="s")
contour(weeks, weeks, tempvarmat, xlab="Days", ylab="Days")
mtext("Temperature Covariance", line=-4, outer=T, cex=1.2)

par(mfrow=c(1,1), pty="s")
persp(weeks, weeks, tempvarmat, theta = 45, phi=45, d=4, shade=0.25, 
      xlab="Days", ylab="Days", zlab="Covariance")
mtext("Temperature Covariance", line=-4, outer=T, cex=1.2)

#  --  plot the temperature correlation function  ---

tempstddev <- sqrt(diag(tempvarmat))
tempcormat <- tempvarmat/outer(tempstddev,tempstddev)
tempcormat <- tempvarmat/(tempstddev %o% tempstddev)

par(mfrow=c(1,1), pty="s")
contour(weeks, weeks, tempcormat, xlab="Days", ylab="Days")
mtext("Temperature Correlation", line=-4, outer=T, cex=1.2)

par(mfrow=c(1,1), pty="s")
persp(weeks, weeks, tempcormat, theta = 45, phi=45, d=4, shade=0.25, 
      xlab="Days", ylab="Days", zlab="Correlation")
mtext("Temperature Correlation", line=-4, outer=T, cex=1.2)

#  --  plot the precipitation variance-covariance bivariate function  ----

precvarbifd <- var.fd(dayprecfd)
precvarmat  <- eval.bifd(weeks,weeks,precvarbifd)

par(mfrow=c(1,1), pty="s")
contour(weeks, weeks, precvarmat, xlab="Days", ylab="Days")
mtext("Precipitation Covariance", line=-4, outer=T, cex=1.2)

par(mfrow=c(1,1), pty="s")
persp(weeks, weeks, precvarmat, theta = 45, phi=45, d=4, shade=0.25,
      xlab="Days", ylab="Days", zlab="Covariance")
mtext("Precipitation Covariance", line=-4, outer=T, cex=1.2)

#  --  plot the precipitation correlation function  ---

precstddev <- sqrt(diag(precvarmat))
preccormat <- precvarmat/(precstddev %o% precstddev)

par(mfrow=c(1,1), pty="s")
contour(weeks, weeks, preccormat, xlab="Days", ylab="Days")
mtext("Precipitation Correlation", line=-4, outer=T, cex=1.2)

par(mfrow=c(1,1), pty="s")
persp(weeks, weeks, preccormat, theta = 45, phi=45, d=4, shade=0.5,
      xlab="Days", ylab="Days", zlab="Covariance")
mtext("Precipitation Correlation", line=-4, outer=T, cex=1.2)

#  -----  compute and plot the covariance between temp. and prec.  --

covbifd <- var.fd(daytempfd, dayprecfd)
covmat  <- eval.bifd(weeks,weeks,covbifd)

par(mfrow=c(1,1), pty="s")
contour(weeks, weeks, covmat, xlab="Weeks", ylab="Weeks")
mtext("Temperature-Precipitation Covariance", line=-4, outer=T, cex=1.2)

par(mfrow=c(1,1), pty="s")
persp(weeks, weeks, covmat, theta = 45, phi=45, d=4, shade=0.25,
      xlab="Weeks", ylab="Weeks", zlab="Covariance")
mtext("Temperature-Precipitation Covariance", line=-4, outer=T, cex=1.2)

#  -----  compute and plot the correlation between temp. and prec.  --

cormat  <- covmat/(tempstddev %o% precstddev)

par(mfrow=c(1,1), pty="s")
contour(weeks, weeks, cormat, xlab="Weeks", ylab="Weeks")
mtext("Temperature-Precipitation Correlation", line=-4, outer=T, cex=1.2)

par(mfrow=c(1,1), pty="s")
persp(weeks, weeks, cormat, theta = 45, phi=45, d=4, shade=0.25, 
      xlab="Weeks", ylab="Weeks", zlab="Correlation")
mtext("Temperature-Precipitation Correlation", line=-4, outer=T, cex=1.2)

#  Suggestion:  redo these plots using this index vector to reorder days
#  so that they run from July 1 to June 30.  In this way, you will see
#  in the middle of the plot that part of the weather year that has the
#  most variation and the most interesting structure.
#  If this reordering is desired, use commands
#  daytime = daytime[JJindex]
#  tempav = tempav[JJindex,]
#  precav = precav[JJindex,]

JJindex = c(182:365, 1:181)



