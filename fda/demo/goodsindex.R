#  -----------------------------------------------------------------------
#                Nondurable Goods Manufacturing Index Analyses
#  -----------------------------------------------------------------------

#  -----------------------------------------------------------------------
#
#                          Overview of the analyses
#
#  Not all functional data come in the form of many replications of 
#  functional observations.  These analyses are of a single very long
#  economic index, with monthly values between 1919 and 2000.  The
#  index exhibits variation on multiple time scales.  
#  Here we see how to use the phase plane plot to look at the seasonal
#  trend in terms of an exchange between two energy states: kinetic 
#  and potential.

#  attach FDA functions

#  Windows ... R

#  Last modified 5 November 2008 by Jim Ramsay
#  Previously 18 August 2007 by Giles Hooker

#  ------------------------------------------------------------------------
#                    Set up the data for analysis
#  ------------------------------------------------------------------------

ndur <- 12
durtime <- rep((0:(ndur-1))/12,81) + rep(1919:1999,1,each=ndur)
goodsrange <- c(1919,2000)
monthlabs <- c("j","f","m","A","M","J","J","A","S","O","N","D")

#  compute log nondurables

lognondur <- log10(nondurables[1:972])

#  compute linear trend

lognondurhat <- lognondur - lsfit(durtime, lognondur)$residuals

#  set up plotting arrangements for one and two panel displays allowing
#  for larger fonts

#  plot the index

par(mfrow=c(1,1), mar=c(5,5,4,2), pty="m")
plot(durtime, nondurables[1:972], type="l", cex=1,
     xlim=goodsrange, ylim=c(0,120),
     xlab="Year", ylab="Nondurable Goods Index")

#  plot the log index and its linear trend

plot(durtime, lognondur, type="l", cex=1,
     xlim=goodsrange, ylim=c(0.7,2.2),
     xlab="Year", ylab="Log10 Nondurable Goods Index" )
lines(durtime, lognondurhat, lty=3)

#  smooth the log data with order 8 splines, knots at data points

lambda     <- 10^(-11)
wtvec      <- rep(1,ndur)
goodsbasis <- create.bspline.basis(goodsrange,81*ndur+4,6)
goodsfdPar <- fdPar(goodsbasis, 4, lambda)

#  smooth the data with smooth.basis.  This takes a fair length of time

lognondursmthfd <- smooth.basis(durtime, lognondur, goodsfdPar)$fd

#  evaluate smooth over a fine mesh of values

durfine <- seq(1919,2000,0.025)
lognondursmthvec <- eval.fd(durfine, lognondursmthfd)
 
#  smooth the data with smooth.Pspline.  This is fast, but
#  does not produce a functional data object for the smooth.

#  smooth.Pspline not available at this time.

#  smoothlist <- smooth.Pspline(durtime, lognondur, w=wtvec, norder=4,
#                             spar=lambda, method=1)
# lognondursmth <- smoothlist$ysmth

#  plot the data and smooth for 1964-1966

index <- durtime >= 1964 & durtime <= 1967
plot(durtime[index], lognondur[index], type="p", cex=1, lwd=2, 
     xlim=c(1964,1967), ylim=c(1.61,1.73),
     xlab="Year", ylab="Log10 NGI",cex.axis=1.5,cex.lab=2)
index <- durfine >= 1964 & durfine <= 1967
lines(durfine[index], lognondursmthvec[index], lwd=2)
lines(c(1965,1965),c(1.61,1.73),lty=2)
lines(c(1966,1966),c(1.61,1.73),lty=2)

#  estimate non-seasonal trend by smoothing with knot
#  at each year.  Use order 4 splines.

longtermbasis <- create.bspline.basis(goodsrange, 83)

longtermfit <- smooth.basis(lognondur, durtime, longtermbasis)$fd

#  compute and plot seasonal trend

seasonfit <- lognondur - eval.fd(durtime, longtermfit)

plot(durtime, seasonfit, type="l", cex=1,
     xlim=goodsrange, ylim=c(-0.08,0.08),
     xlab="Year", ylab="Seasonal Trend")
lines(goodsrange, c(0,0), lty=2)

#  plot a sine to illustrate kinetic and potential energy

tval <- seq(0,1,0.01)
xval <- sin(2*pi*tval)

plot(tval, xval, type="l", cex=1,
	  xlab="Time", ylab="Horizontal Position")
lines(c(0,1), c(0,0), lty=2)

#  a phase-plane plot of the sine

Dxval  <-  (2*pi)*cos(2*pi*tval)
D2xval <- -(2*pi)^2*sin(2*pi*tval)

plot(Dxval, D2xval, type="l", 
     xlab="Velocity", ylab="Acceleration")
lines(c(-2*pi,2*pi), c(0,0), lty=2) 
lines(c(0,0), c(-(2*pi)^2,(2*pi)^2), lty=2)
 
#  ----------  phase-plane plots of the goods index  -------------------

#  -----------------  1964  -----------------------

nyr       <-  1   #  number of startyr to plot
startyr   <- 64   #   starting year
yearindex <- startyr + 1900 - 1919  # indices of year

#  This is the code that carries out the phase-plane plotting.

index    <- (1:(nyr*12+1)) + yearindex*12
durrange <- c(min(durtime[index]), max(durtime[index]))
durfine  <- seq(durrange[1],durrange[2],len=401)
durcrse  <- seq(durrange[1],durrange[2],len=nyr*12+1)
D1f <- eval.fd(durfine, lognondursmthfd, 1)
D2f <- eval.fd(durfine, lognondursmthfd, 2)
D1c <- eval.fd(durcrse, lognondursmthfd, 1)
D2c <- eval.fd(durcrse, lognondursmthfd, 2)

plot(c(-.75,.75),c(-12,12), type="n", cex=1,
     xlab="Velocity",  ylab="Acceleration")
lines(c(0,0), c(-12,12), lty=2) 
lines(c(-.75,.75), c(0,0), lty=2)
x0 <- durrange[1]
indj <- durfine >= x0 & durfine <= x0+1
lines(D1f[indj], D2f[indj])
indexc <-  1:12 
for (k in 1:12) {
    points(D1c[indexc[k]], D2c[indexc[k]])
    text(D1c[indexc[k]]+0.05, D2c[indexc[k]]+0.5, monthlabs[k])
}
title(paste("Year", startyr + 1900))


