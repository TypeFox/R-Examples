#  -----------------------------------------------------------------------
#                      Melanoma Incidence data
#  -----------------------------------------------------------------------
#
#                          Overview of the analyses
#
#  The data are numbers of melanoma cases per 100,000 in the state of 
#  Connecticut for the years 1936 to 1972.  The data show two kinds
#  of trend:  a nearly linear long-term trend over this period,
#  and a cycle of around 10 years corresponding to the
#  sunspot cycle.  
#  These analyses are designed to contrast two types of
#  smoothing techniques.
#  The first and simplest is aa straightforward
#  smooth of the data using B-spline basis functions defined
#  by knots at each data point and using a roughness penalty
#  that penalizes the total size of the second derivative or
#  curvature.
#  The second smooth defines a fourth order linear differential
#  operator that penalizes departure from a linear trend plus
#  sinusoidal trend defined by a period of 9.67 years.  This 
#  more sophisticated penalty highlights departures from this
#  type of baseline trend.  
#  -----------------------------------------------------------------------

#  Last modified  21 March 2006

#  attach the FDA functions
#  input the data  

tempmat <- t(matrix(scan("../data/melanoma.txt", 0), 3, 37))

year  <- tempmat[,2]
mela  <- tempmat[,3]
nyear <- length(year)
rng   <- c(min(year),max(year))

melahat1 <- mela - lsfit(year, mela)$residual
sse1 <- sum((mela-melahat1)^2)

xmat2 <- cbind(year,sin(0.65*year))
melahat2 <- mela - lsfit(xmat2, mela)$residual
sse2 <- sum((mela-melahat2)^2)

par(mfrow=c(1,1), mar=c(5,5,4,2), pty="m")
plot (year, mela, type="p", cex=1,
      xlab="Year", ylab="Cases per 100,000")
lines(year, melahat1, lty=4)
lines(year, melahat2, lty=1)

lnmela <- log10(mela)

lnmelahat1 <- lnmela - lsfit(year, lnmela)$residual
sse1 <- sum((lnmela-lnmelahat1)^2)

lnmelahat2 <- lnmela - lsfit(xmat2, lnmela)$residual
sse2 <- sum((lnmela-lnmelahat2)^2)

plot (year, lnmela, type="p", cex=1, ylim=c(-0.1,0.8),
      xlab="Year", ylab="Log_10 Cases per 100,000")
lines(year, lnmelahat1, lty=4)
lines(year, lnmelahat2, lty=1)

#  -----------------------------------------------------------------------
#             smooth data using B-splines
#  -----------------------------------------------------------------------

#  set up the basis with a knot at every year

knots  <- year
nbasis <- nyear + 4
norder <- 6
basisobj  <- create.bspline.basis(rng, nbasis, norder, knots)

#  smooth the data by penalizing the second derivative

Lfdobj <- 2
lambda <- 1
melafdPar <- fdPar(basisobj, Lfdobj, lambda)

melafd <- smooth.basis(year, mela, melafdPar)$fd

#  plot the data and the smooth

plotfit.fd(mela, year, melafd)

#  plot the residuals

plotfit.fd(mela, year, melafd, residual=TRUE)

#  set up operator to remove sinusoid plus trend

omega   <- 0.65
Lfdobj  <- vec2Lfd(c(0, 0, omega^2, 0), rng)
lambda  <- 1e-2
melafdPar <- fdPar(basisobj, Lfdobj, lambda)

#  smooth the data

melafd <- smooth.basis(year, mela, melafdPar)$fd

#  plot the results

plotfit.fd(mela, year, melafd)







