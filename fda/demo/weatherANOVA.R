#  Analysis variance for the temperature and precipitation data with
#  climate zone being the grouping factor

#  This file is intended to be used after the commands in files
#  weathersetup.R and weathersmooth.R have been executed.

#  Many other interesting ways of describing the data and plotting results
#  can be found in the file canadian-weather.R, set up in 2008 by 
#  Spencer Graves.

#  Last modified 17 November 2008

load("weatherdata")

station <- weatherdata$station

#  names for (climate zones

zonenames <- c("Canada  ", "Atlantic", "Pacific ", "Contintl", "Arctic  ")

#  indices for (weather stations in each of four climate zones

#  set up indices that order the stations from east to west to north

atlindex <-  1:15
conindex <- 16:24
pacindex <- 25:31
artindex <- 32:35

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

#  ---------------------------------------------------------------
#             Predicting temperature from climate zone
#  ---------------------------------------------------------------

#  We use a more compact basis here to avoid storage allocation
#  problems that arise when 365 basis functions are used.

dayrange   <- c(0,365)
daybasis65 <- create.fourier.basis(dayrange, 65)
daytime    <- weatherdata$daytime
tempav     <- weatherdata$tempav
smoothList <- smooth.basis(daytime, tempav, daybasis65)
daytempfd  <- smoothList$fd
tempy2cMap <- smoothList$y2cMap

daytempfd$fdnames <- list(NULL, station, NULL)

#  revise YFDOBJ by adding a zero function

coef   <- daytempfd$coefs
coef36 <- cbind(coef,matrix(0,65,1))
daytempfd$coefs <- coef36

p <- 5
xfdlist <- vector("list",p)
for (j in 1:p) xfdlist[[j]] <- zmat[,j]

#  set up the basis for (the regression functions

nbetabasis <- 13
betabasis  <- create.fourier.basis(dayrange, nbetabasis)

#  set up the harmonic acceleration operator

harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), dayrange)

#  set up the functional parameter object for (the regression fns.

betafd    <- fd(matrix(0,nbetabasis,1), betabasis)
estimate  <- T
lambda    <- 0
betafdPar <- fdPar(betafd, harmaccelLfd, lambda, estimate)

betalist <- vector("list",p)
for (j in 1:p) betalist[[j]] <- betafdPar

#  compute regression coefficient functions and
#  predicted functions

fRegressList <- fRegress(daytempfd, xfdlist, betalist)

# try to get a cross-validated score

res.cv = fRegress.CV(daytempfd,xfdlist,betalist,whichobs=1:35)


#  plot regression functions

betaestlist <- fRegressList$betaestlist

par(mfrow=c(3,2))
for (j in 1:p) {
	betaestParfdj <- betaestlist[[j]]
	plot(betaestParfdj$fd, xlab="Day", ylab="Temp.",
	     main=zonenames[j])
}

#  set up predicted functions

yhatfdobj <- fRegressList$yhatfdobj

#  compute residual matrix and get covariance of residuals

yhatmat    <- eval.fd(daytime, yhatfdobj)
ymat       <- eval.fd(daytime, daytempfd)
tempresmat <- ymat[,1:35] - yhatmat[,1:35]
SigmaE     <- var(t(tempresmat))

#  plot covariance surface for errors

par(mfrow=c(1,1))
contour(SigmaE, xlab="Day", ylab="Day")
lines(dayrange,dayrange,lty=4)

#  plot standard deviation of errors

par(mfrow=c(1,1), mar=c(5,5,3,2), pty="m", ask=F)
stddevE <- sqrt(diag(SigmaE))
plot(daytime, stddevE, type="l",
     xlab="Day", ylab="Standard error (deg C)")

#  Repeat regression, this time outputting results for
#  confidence intervals

stderrList <- fRegress.stderr(fRegressList, tempy2cMap, SigmaE)

betastderrlist <- stderrList$betastderrlist

#  plot regression function standard errors

par(mfrow=c(2,3))
for (j in 1:p) {
	betastderrj <- eval.fd(daytime, betastderrlist[[j]])
	plot(daytime, betastderrj,
	        type="l",lty=1, xlab="Day", ylab="Reg. Coeff.",
	        main=zonenames[j])
}

#  plot regression functions with confidence limits

par(mfrow=c(2,3))
for (j in 1:p) {
	betafdParj  <- betaestlist[[j]]
	betafdj     <- betafdParj$fd
	betaj       <- eval.fd(daytime, betafdj)
	betastderrj <- eval.fd(daytime, betastderrlist[[j]])
	matplot(daytime, cbind(betaj, betaj+2*betastderrj, betaj-2*betastderrj),
	        type="l",lty=c(1,4,4), xlab="Day", ylab="Reg. Coeff.",
	        main=zonenames[j])
}

#  set up a functional data object for the temperature residuals
 
lambda   <- 1e5
fdParobj <- fdPar(daybasis65, harmaccelLfd, lambda)
 
smoothList <- smooth.basis(daytime, tempresmat, fdParobj)

tempresfdobj <- smoothList$fd

#  plot temperature residuals

par(mfrow=c(1,1))
plot(tempresfdobj, ask=F)
 
#  ---------------------------------------------------------------
#             Predicting precipitation from climate zone
#  ---------------------------------------------------------------

#  We use a more compact basis here to avoid storage allocation
#  problems that arise when 365 basis functions are used.

daybasis21 <- create.fourier.basis(dayrange, 21)
daytime    <- weatherdata$daytime
precav     <- weatherdata$precav
smoothList <- smooth.basis(daytime, precav, daybasis21)
dayprecfd  <- smoothList$fd
precy2cMap <- smoothList$y2cMap

#  revise YFDOBJ by adding a zero function

coef   <- dayprecfd$coefs
coef36 <- cbind(coef,matrix(0,21,1))
dayprecfd$coefs <- coef36

p <- 5
xfdlist <- vector("list",p)
for (j in 1:p) xfdlist[[j]] <- zmat[,j]

#  set up the basis for (the regression functions

nbetabasis <- 13
betabasis  <- create.fourier.basis(dayrange, nbetabasis)

#  set up the functional parameter object for (the regression fns.

betafd    <- fd(matrix(0,nbetabasis,1), betabasis)
estimate  <- T
lambda    <- 0
betafdPar <- fdPar(betafd, harmaccelLfd, lambda, estimate)

betalist <- vector("list",p)
for (j in 1:p) betalist[[j]] <- betafdPar

#  compute regression coefficient functions and
#  predicted functions

fRegressList <- fRegress(dayprecfd, xfdlist, betalist)

#  plot regression functions

betaestlist <- fRegressList$betaestlist

par(mfrow=c(3,2))
for (j in 1:p) {
	betaestParfdj <- betaestlist[[j]]
	plot(betaestParfdj$fd, xlab="Day", ylab="Prec.",
	     main=zonenames[j])
}

#  set up predicted functions

yhatfdobj <- fRegressList$yhatfdobj

#  compute residual matrix and get covariance of residuals

yhatmat    <- eval.fd(daytime, yhatfdobj)
ymat       <- eval.fd(daytime, dayprecfd)
precresmat <- ymat[,1:35] - yhatmat[,1:35]
SigmaE     <- var(t(precresmat))

#  plot covariance surface for errors

par(mfrow=c(1,1))
contour(SigmaE, xlab="Day", ylab="Day")
lines(dayrange,dayrange,lty=4)

#  plot standard deviation of errors

par(mfrow=c(1,1), mar=c(5,5,3,2), pty="m", ask=F)
stddevE <- sqrt(diag(SigmaE))
plot(daytime, stddevE, type="l",
     xlab="Day", ylab="Standard error (mm)")

#  Repeat regression, this time outputting results for
#  confidence intervals

stderrList <- fRegress.stderr(fRegressList, precy2cMap, SigmaE)

betastderrlist <- stderrList$betastderrlist

#  plot regression function standard errors

par(mfrow=c(2,3))
for (j in 1:p) {
	betastderrj <- eval.fd(daytime, betastderrlist[[j]])
	plot(daytime, betastderrj,
	        type="l",lty=1, xlab="Day", ylab="Reg. Coeff.",
	        main=zonenames[j])
}

#  plot regression functions with confidence limits

par(mfrow=c(2,3))
for (j in 1:p) {
	betafdParj  <- betaestlist[[j]]
	betafdj     <- betafdParj$fd
	betaj       <- eval.fd(daytime, betafdj)
	betastderrj <- eval.fd(daytime, betastderrlist[[j]])
	matplot(daytime, cbind(betaj, betaj+2*betastderrj, betaj-2*betastderrj),
	        type="l",lty=c(1,4,4), xlab="Day", ylab="Reg. Coeff.",
	        main=zonenames[j])
}

#  set up a functional data object for the precipitation residuals
 
lambda   <- 1e5
fdParobj <- fdPar(daybasis21, harmaccelLfd, lambda)
 
smoothList <- smooth.basis(daytime, precresmat, fdParobj)

precresfdobj <- smoothList$fd

#  plot precipitation residuals

par(mfrow=c(1,1))
plot(precresfdobj, ask=F)
 
#  save temperature and precipitation residual objects

weatherANOVA <- list(tempresfd=tempresfdobj, precpresfd=precresfdobj)

save(weatherANOVA, file="weatherANOVA")



