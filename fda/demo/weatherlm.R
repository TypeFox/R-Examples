#  Some linear models for the daily weather data

#  This file is intended to be used after the commands in files
#  weathersetup.R and weathersmooth.R have been executed.

#  Many other interesting ways of describing the data and plotting results
#  can be found in the file canadian-weather.R, set up in 2008 by 
#  Spencer Graves.

#  Last modified 17 November 2008

load("weatherfd")
load("weatherdata")
load("weatherANOVA")

#  ---------  create fd objects for temp. and prec. ---------------

tempfdPar <- weatherfd$tempfdPar
precfdPar <- weatherfd$precfdPar

daytempfd <- tempfdPar$fd
dayprecfd <- precfdPar$fd

harmaccelLfd <- tempfdPar$Lfd

station <- weatherdata$station

#  We use a more compact basis here to avoid storage allocation
#  problems that arise when 365 basis functions are used.

daybasis65 <- create.fourier.basis(c(0,365), 65)
daytime    <- weatherdata$daytime
precav     <- weatherdata$precav
smoothList <- smooth.basis(daytime, precav, daybasis65)
dayprecfd  <- smoothList$fd
precy2cMap <- smoothList$y2cMap

#  -----------------------------------------------------------------------
#         predict log precipitation from climate zone and temperature
#  -----------------------------------------------------------------------

#  Be sure to run previous analysis predicting temperature from
#  climate zone before running this example.

#  set up functional data object for log precipitation

logprecmat <- log10(eval.fd(daytime, dayprecfd))

lndayprecfd <- smooth.basis(daytime, logprecmat, daybasis65)$fd

lndayprecfd$fdnames <- list(Day = NULL, Station = station, Log10DegC = NULL)

#  plot log precipitation functions

par(mfrow=c(1,1), pty="m")
plot(lndayprecfd)
title("Log Precipitation Functions")

#  revise LOGPREDFD by adding a zero function

coef   <- lndayprecfd$coefs
nbasis <- dim(coef)[1]
coef36 <- cbind(coef,matrix(0,nbasis,1))
lndayprecfd$coefs <- coef36

#  set up the XFDLIST list

p <- 6
xfdlist <- vector("list",p)

#  names for (climate zones

zonenames <- c("Canada  ",
               "Atlantic", "Pacific ", "Contintl", "Arctic  ")

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

#  load first five members with columns of design matrix

for (j in 1:5) xfdlist[[j]] <- matrix(zmat[,j],1,36)

#  extend temperature residual functions to include
#  zero function

tempresfd <- weatherANOVA$tempresfd

coef   <- tempresfd$coefs
nbasis <- dim(coef)[1]
coef36 <- cbind(coef,matrix(0,nbasis,1))
tempresfd$coefs <- coef36

#  add TEMPRFDOBJ to the set of predictors

xfdlist[[6]]  <- tempresfd
betalist[[6]] <- betafdPar

#  set up the basis for (the regression functions

nbetabasis <- 13
betabasis  <- create.fourier.basis(dayrange, nbetabasis)

#  set up the functional parameter object for (the regression fns.

betafd    <- fd(matrix(0,nbetabasis,p), betabasis)
estimate  <- T
lambda    <- 0
betafdPar <- fdPar(betafd, harmaccelLfd, lambda, estimate)
for (j in 1:p) betalist[[j]] <- betafdPar

#  compute regression coefficient functions and
#  predicted functions

fRegressList <- fRegress(lndayprecfd, xfdlist, betalist)

betaestlist <- fRegressList$betaestlist
yhatfdobj   <- fRegressList$yhatfdobj

#  plot regression functions

prednames <- c(zonenames, "tempres ")
par(mfrow=c(2,3),pty="s")
for (j in 1:p) {
    betaParfdj <- betaestlist[[j]]
    betafdj    <- betaParfdj$fd
    plot(betafdj)
    title(prednames[j])
}

#  plot predicted functions

par(mfrow=c(1,1), pty="m")
plot(yhatfdobj)

#  compute residual matrix and get covariance of residuals

yhatmat    <- eval.fd(daytime, yhatfdobj)
ymat       <- eval.fd(daytime, lndayprecfd)
lnprecrmat <- ymat[,1:35] - yhatmat[,1:35]
SigmaE     <- var(t(lnprecrmat))

dayrange <- c(0,365)
par(mfrow=c(1,1))
contour(SigmaE, xlab="Day", ylab="Day")
lines(dayrange,dayrange,lty=4)

#  repeat regression analysis to get confidence intervals

stderrList <- fRegress.stderr(fRegressList, precy2cMap, SigmaE)

betastderrlist <- stderrList$betastderrlist

#  plot regression function standard errors

par(mfrow=c(2,3), pty="s")
for (j in 1:p) {
	betastderrfdj <- betastderrlist[[j]]
	betastderrj <- eval.fd(daytime, betastderrfdj)
	plot(daytime, betastderrj,
	        type="l",lty=1, xlab="Day", ylab="Reg. Coeff.",
	        main=prednames[j])
}

#  plot regression functions with confidence limits

par(mfrow=c(2,3), pty="s")
for (j in 1:p) {
	betafdParj  <- betaestlist[[j]]
	betafdj     <- betafdParj$fd
	betaj       <- eval.fd(daytime, betafdj)
	betastderrfdj <- betastderrlist[[j]]
	betastderrj   <- eval.fd(daytime, betastderrfdj)
	matplot(daytime, cbind(betaj, betaj+2*betastderrj, betaj-2*betastderrj),
	        type="l",lty=c(1,4,4), xlab="Day", ylab="Reg. Coeff.",
	        main=prednames[j])
}

#  ---------------------------------------------------------------
#      log annual precipitation predicted by temperature profile
#  ---------------------------------------------------------------

#  set up log10 total precipitation

annualprec <- log10(apply(precav,2,sum))

#  We use a more compact basis here to avoid storage allocation
#  problems that arise when 365 basis functions are used.

daybasis65 <- create.fourier.basis(c(0,365), 65)
daytime    <- weatherdata$daytime
tempav     <- weatherdata$tempav
smoothList <- smooth.basis(daytime, tempav, daybasis65)
daytempfd  <- smoothList$fd
tempy2cMap <- smoothList$y2cMap

#  set up the covariates, the first the constant, and the second
#  temperature

p <- 2
constantfd <- fd(matrix(1,1,35), create.constant.basis(dayrange))

xfdlist <- vector("list",2)
xfdlist[[1]] <- constantfd
xfdlist[[2]] <- daytempfd

#  set up the functional parameter object for (the regression fns.
#  the smoothing parameter for (the temperature function
#  is obviously too small here, and will be revised below by
#  using cross-validation.

betalist   <- vector("list",2)
#  set up the first regression function as a constant
betabasis1 <- create.constant.basis(dayrange)
betafd1    <- fd(0, betabasis1)
betafdPar1 <- fdPar(betafd1)
betalist[[1]] <- betafdPar1
#  set up the second with same basis as for (temperature
#  35 basis functions would permit a perfect fit to the data
nbetabasis  <- 35
betabasis2  <- create.fourier.basis(dayrange, nbetabasis)
betafd2     <- fd(matrix(0,nbetabasis,1), betabasis2)
lambda      <- 10
betafdPar2  <- fdPar(betafd2, harmaccelLfd, lambda)
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
#    SSE.CV[ilam]   <- fRegress.CV(annualprec, xfdlist, betalisti)$SSE.CV
    SSE.CV[ilam] <- fRegress(annualprec, xfdlist, betalisti)$OCV
    print(c(ilam, loglam[ilam], SSE.CV[ilam]))
}

plot(loglam, SSE.CV, type="b",
	  xlab="log_{10} smoothing parameter lambda",
	  ylab="Cross-validation score")

#  analysis with minimum CV smoothing

lambda        <- 10^12.5
betafdPar2    <- fdPar(betafd2, harmaccelLfd, lambda)
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
plot(betafd, xlab="Day", ylab="Regression function")
title("Regression coefficient for temperature")

#  plot the fit

par(mfrow=c(1,1), pty="m")
plot (annualprechat, annualprec, type="p")
lines(annualprechat, annualprechat, lty=2)

#  compute squared multiple correlation

covmat <- var(cbind(annualprec, annualprechat))
Rsqrd <- covmat[1,2]^2/(covmat[1,1]*covmat[2,2])
#   0.75

#  compute SigmaE

resid  <- annualprec - annualprechat
SigmaE <- mean(resid^2)
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
betavec         <- eval.fd(daytime, betafd)
betastderrvec   <- eval.fd(daytime, betastderrfd)

betaplotmat <- cbind(betavec, betavec+2*betastderrvec,
                              betavec-2*betastderrvec)

matplot(daytime, betaplotmat, type="l", lty=c(1,4,4),
        xlab="Day", ylab="Temperature Reg. Coeff.")
lines(dayrange,c(0,0),lty=2)

#  ---------------------------------------------------------------
#         predict log precipitation from temperature
#  ---------------------------------------------------------------

#  set up a smaller basis using only 65 Fourier basis functions
#  to save some computation time

smallnbasis <- 65
smallbasis  <- create.fourier.basis(dayrange, smallnbasis)

tempfd      <- smooth.basis(daytime, tempav, smallbasis)$fd

#  change 0's to 0.05 mm in precipitation data

prectmp <- precav
for (j in 1:35) {
    index <- prectmp[,j] == 0
    prectmp[index,j] <- 0.05
}

#  work with log base 10 precipitation

logprecmat <- log10(prectmp)

#  set up functional data object for log precipitation

lndayprecfd <- smooth.basis(daytime, logprecmat, smallbasis)$fd
lndayprecfdnames <- vector("list",3)
lndayprecfdnames[[1]] <- "Days"
lndayprecfdnames[[2]] <- "Station"
lndayprecfdnames[[3]] <- "log.{10} mm"
lndayprecfd$fdnames <- lndayprecfdnames

#  plot precipitation functions

par(mfrow=c(1,1), pty="m")
plot(lndayprecfd)
title("Log Precipitation Functions")

#  set up smoothing levels for (s (xLfd) and for (t (yLfd)

xLfdobj <- harmaccelLfd
yLfdobj <- harmaccelLfd
xlambda <- 1e9
ylambda <- 1e7

#  compute the linear model

wtvec <- matrix(1,35,1)

linmodstr <- linmod(daytempfd, lndayprecfd, wtvec,
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

lnprecmeanfd <- mean(lndayprecfd)
lnprechat0   <- eval.fd(weeks, lnprecmeanfd)

#  Plot actual observed, fitted, and mean log precipitation for
#      each weather station,

lnprecmat    <- eval.fd(weeks, lndayprecfd)
lnprechatmat <- eval.fd(weeks, lnprechatfd)
plotrange    <- c(min(lnprecmat),max(lnprecmat))
par(ask=T)
for (i in 1:35) {
    lnpreci    <- eval.fd(lndayprecfd[i],    weeks)
    lnprechati <- eval.fd(lnprechatfd[i], weeks)
    SSE <- sum((lnpreci-lnprechati)^2)
    SSY <- sum((lnpreci-lnprechat0)^2)
    RSQ <- (SSY-SSE)/SSY
    plot(weeks, lnprecmat[,i], type="l", lty=4, ylim=plotrange,
         xlab="Day", ylab="Log Precipitation",
         main=paste(place[i],"  R^2 =",round(RSQ,3)))
    lines(weeks, lnprechatmat[,i])
}
par(ask=F)

