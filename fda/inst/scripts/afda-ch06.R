###
###
### Ramsey & Silverman (2002) Applied Functional Data Analysis
###                           (Springer)
###
### ch. 6.  Human growth 
###
library(fda)

##
## sec. 6.1.  Introduction
##

##
## sec. 6.2.  Height measurements at three scales
##
str(growth)

# pp.  84-85.  Figure 6.1.  the first 10 females
# of the Berkeley growth study
op <- par(cex=1.1)
with(growth, matplot(age, hgtf[, 1:10], type="b", pch="o", 
                     ylab="Height (cm.)") )
par(op)
# Figure 6.2.  Heights of one boy during one school year

# A monotone smooth requires some effort.  
#  (1) set up the basis

nbasis.onechild   <- 33
# Establish a B-spline basis
# with nbasis.onechild basis function
hgtbasis <- with(onechild, 
  create.bspline.basis(range(day), nbasis.onechild))

tst <- create.bspline.basis(onechild$day)
#  set up the functional data object for W <- log Dh

# Start by creating a functional 0 from hgtbasis
cvec0 <- rep(0,nbasis.onechild)
Wfd0   <- fd(cvec0, hgtbasis)

#  set parameters for the monotone smooth
# with smoothing lambda = 1e-1 or 1e-12 

WfdPar.1  <- fdPar(Wfd0, 2, lambda=.1)
WfdPar.12  <- fdPar(Wfd0, 2, lambda=1e-12)

#  --------------   carry out the monotone smooth  ---------------

# The monotone smooth is
# beta[1]+beta[2]*integral(exp(Wfdobj)),
#     where Wfdobj is a functional data object
smoothList.1   <- with(onechild, 
  smooth.monotone(x=day, y=height, WfdParobj=WfdPar.1) )
smoothList.12   <- with(onechild, 
  smooth.monotone(x=day, y=height, WfdParobj=WfdPar.12) )

#str(smoothList)
#attach(smoothList)

# Create a fine grid at which to evaluate the smooth
dayfine  <- with(onechild, seq(day[1],day[length(day)],len=151))

# eval.monfd = integral(exp(Wfdobj))
# This is monotonically increasing, since exp(Wfdobj)>0.  
#yhat     <- with(smoothList,
#             beta[1] + beta[2]*eval.monfd(onechild$day, Wfdobj))
yhatfine.1 <- with(smoothList.1, 
       beta[1] + beta[2]*eval.monfd(dayfine, Wfdobj))
yhatfine.12 <- with(smoothList.12, 
       beta[1] + beta[2]*eval.monfd(dayfine, Wfdobj))

plot(onechild, ylab="Height (cm.)") # raw data
lines(dayfine, yhatfine.1, lwd=2) # lambda=0.1:  reasonable
lines(dayfine, yhatfine.12, lty=2)
# lambda=1e-12:too close to a straight line

# p.  86, Figure 6.3.  Growth of the length of the tibia of a newborn;  
# data not available





##
## sec. 6.3.  Velocity and acceleration
##

# p. 87, Figure 6.4.  Estimated growth rate of the first girl 

nage <- length(growth$age)
norder.growth <- 6
nbasis.growth <- nage + norder.growth - 2
# 35 
rng.growth <- range(growth$age)
# 1 18 
wbasis.growth <- create.bspline.basis(rng.growth,
                   nbasis.growth, norder.growth,
                   growth$age)

#  starting values for coefficient

cvec0.growth <- matrix(0,nbasis.growth,1)
Wfd0.growth  <- fd(cvec0.growth, wbasis.growth)

Lfdobj.growth    <- 3          #  penalize curvature of acceleration
lambda.growth    <- 10^(-0.5)  #  smoothing parameter
growfdPar <- fdPar(Wfd0.growth, Lfdobj.growth, lambda.growth)

#  ---------------------  Now smooth the data  --------------------

smoothGirl1 <- with(growth, smooth.monotone(x=age,
        y=hgtf[, 1], WfdParobj=growfdPar, conv=0.001, active=TRUE, 
        dbglev=0) )
#*** The default active = c(FALSE, rep(TRUE, ncvec-1))
#    This tells smooth.monotone to estimate
#    force the first coeffeicient of W to 0,
#    and estimate only the others.
#*** That starts velocity unrealistically low.
#*** Fix this with active=TRUE 
#    to estimate all elements of W.  

# Create a fine grid at which to evaluate the smooth
agefine  <- with(growth, seq(age[1], age[nage], len=151))

# Lfdobj = 1 for first derivative, growth rate or velocity 
smoothG1.1 <- with(smoothGirl1,
       beta[2]*eval.monfd(agefine, Wfdobj, Lfdobj=1))

plot(agefine, smoothG1.1, type="l",
     xlab="Year", ylab="Growth velocity (cm/year)")
axis(3, labels=FALSE)
axis(4, labels=FALSE)


# p. 88, Figure 6.5, Estimated growth velocity of a 10-year old boy
str(onechild)

nDays <- dim(onechild)[1]
norder.oneCh <- 6
nbasis.oneCh <- nDays+norder.oneCh-2
rng.days <- range(onechild$day)

# B-spline basis 
wbasis.oneCh <- create.bspline.basis(rng.days,
       nbasis.oneCh, norder.oneCh, onechild$day)

# starting values for coefficients
cvec0.oneCh <- matrix(0, nbasis.oneCh, 1)
Wfd0.oneCh <- fd(cvec0.oneCh, wbasis.oneCh)

Lfdobj.oneCh    <- 3
#  penalize curvature of acceleration
growfdPar.oneCh100 <- fdPar(Wfd0.oneCh, Lfdobj.oneCh,
                 lambda=100 )

# now smooth the data

zmat.oneCh <- matrix(1, nDays, 1)
smoothOneCh100 <- with(onechild, smooth.monotone(x=day,
       y=height, WfdParobj=growfdPar.oneCh100, 
       conv=0.001, active=TRUE, dbglev=0) )

dayFine <- with(onechild, seq(day[1], day[nDays], len=151))

sm.OneCh100 <- with(smoothOneCh100,
   beta[2]*eval.monfd(dayFine, Wfdobj, Lfdobj=1) )

plot(dayFine, sm.OneCh100, type="l", xlab="day",
     ylab="Growth velocity (cm/day)" )
axis(3, labels=FALSE)
axis(4, labels=FALSE)
# Close to Figure 6.5
# closer than with lambda = 10 or 1000 

# p. 88, Figure 6.6.  Estimated growth velocity of a baby

# Data not available.






# p. 89, Figure 6.7.
#Estimated growth acceleration for 10 girls in the Berkeley Growth Study
# Similar to Figure 6.4, but acceleration not velocity
# and for 10 girls not one

nage <- length(growth$age)
norder.growth <- 6
nbasis.growth <- nage + norder.growth - 2
# 35 
rng.growth <- range(growth$age)
# 1 18 
wbasis.growth <- create.bspline.basis(rng.growth,
                   nbasis.growth, norder.growth,
                   growth$age)
str(wbasis.growth)
#  starting values for coefficient

cvec0.growth <- matrix(0,nbasis.growth,1)
Wfd0.growth  <- fd(cvec0.growth, wbasis.growth)

Lfdobj.growth    <- 3          #  penalize curvature of acceleration

#  ---------------------  Now smooth the data  --------------------

# Create a fine grid at which to evaluate the smooth
nptsFine <- 151

#############################################
# Experimented with other values of lambda.
# Relative to Figure 6.7,
# lambda = 0.01 undersmoothed, while 0.1 oversmoothed.
# NOTE:  This script was prepared years after the original
# figures were prepared, without knowledge of
# exactly how the original figures were prepared.

lambda.gr2.3 <- .03
growfdPar2.3 <- fdPar(Wfd0.growth, Lfdobj.growth, lambda.gr2.3)

ncasef <- 10
smoothGirls2.3 <- vector("list", ncasef)
for(icase in 1:ncasef){
  smoothGirls2.3[[icase]] <- with(growth, smooth.monotone(x=age,
      y=hgtf[, icase], WfdParobj=growfdPar2.3, conv=0.001,
      active=TRUE, dbglev=0) )
  cat(icase, "")
}
    
agefine  <- with(growth, seq(age[1], age[nage], len=nptsFine))
smoothGirlsAcc2.3 <- array(NA, dim=c(nptsFine, ncasef))
for(icase in 1:ncasef)
  smoothGirlsAcc2.3[, icase] <- with(smoothGirls2.3[[icase]],
     beta[2]*eval.monfd(agefine, Wfdobj, Lfdobj=2) )

op <- par(cex=1.5)
matplot(agefine, smoothGirlsAcc2.3, type="l", ylim=c(-4, 2),
        xlab="age", ylab="Growth acceleration (cm/year^2)",
        bty="n")
lines(agefine, apply(smoothGirlsAcc2.3, 1, mean), lwd=3)
abline(h=0, lty="dashed")
par(op)
# Good match

#############################################

##
## sec. 6.4.  An equation for growth 
##

# p. 91, Figure 6.8.  Relative acceleration and its integral

#Estimated growth acceleration for 10 girls in the Berkeley Growth Study
# Similar to Figure 6.7, but acceleration / velocity

smoothGirlsVel2.3 <- array(NA, dim=c(nptsFine, ncasef))
for(icase in 1:ncasef)
  smoothGirlsVel2.3[, icase] <- with(smoothGirls2.3[[icase]],
     beta[2]*eval.monfd(agefine, Wfdobj, Lfdobj=1) )

smoothGirlsRelAcc2.3 <- (smoothGirlsAcc2.3 /
                         smoothGirlsVel2.3) 

matplot(agefine, smoothGirlsRelAcc2.3, type="l", ylim=c(-2, 0.5),
        xlab="Year", ylab="w(t)")
abline(h=0, lty="dashed")

diff(range(quantile(diff(agefine))))
d.agefine <- mean(diff(agefine))
smoothGirlsIntRelAcc2.3 <- (apply(smoothGirlsRelAcc2.3, 2, cumsum)
                            * d.agefine)
str(smoothGirlsIntRelAcc2.3)

matplot(agefine, smoothGirlsIntRelAcc2.3, type="l", ylim=c(-4, 0.5),
        xlab="Year", ylab="W(t)")
abline(h=0, lty="dashed")

##
## sec. 6.5.  Timing or phase variation in growth 
##

# p. 92, Figure 6.9.  Acceleration curves
#  differing in amplitude only and phase only

# This figure was conceptual, not real data.

# Recreate roughly by reading numbers off the figure
Fig6.9 <- cbind(Age=3:20,
                accel=c(-0.6, -0.5, -0.41, -0.4, -0.4, -0.405,
                  -0.4, -0.25, 0.2, 1.9, 2.8, -1.8, -4.7,
                  -3.2, -0.7, -0.2, -0.02, 0) ) 
Fig6.9basis6 <- create.bspline.basis(Fig6.9[, 1], norder=5)

# If Lfdobj = 0 (default), smooth.basisPar goes crazy
# between points
# Lfdobj = 3 works OK with lambda = 0.01
Fig6.9s3.01 <- smooth.basisPar(argvals=Fig6.9[, 1], y=Fig6.9[, 2],
               fdobj=Fig6.9basis6, 3, .01)
plot(Fig6.9s3.01$fd, ylim=range(Fig6.9[, 2]))

Fig6.9age <- seq(5, 20, length=101)
Fig6.9a <- eval.fd(Fig6.9age, Fig6.9s3.01$fd)
plot(Fig6.9age, 1.2*Fig6.9a, type="l", xlab="Age", ylab="")
lines(Fig6.9age, Fig6.9a, lty=2, col=2)
lines(Fig6.9age, 0.8*Fig6.9a, lty=3, col=3)
title("Amplitude variation")

Fig6.9b1 <- seq(3, 18, length=101)
Fig6.9b2 <- seq(4, 19, length=101)
Fig6.9b1. <- eval.fd(Fig6.9b1, Fig6.9s3.01$fd)
Fig6.9b2. <- eval.fd(Fig6.9b2, Fig6.9s3.01$fd)

plot(Fig6.9age, Fig6.9a, type="l", xlab="Age", ylab="")
lines(Fig6.9age, Fig6.9b1., lty=2, col=2)
lines(Fig6.9age, Fig6.9b2., lty=3, col=3)
title("Phase variation")

# p. 93, Figure 6.10.  Time warping functions
#  for 10 Berkeley girls

#  ---------------------------------------------------------------------
#            Register the velocity curves for the girls
#  ---------------------------------------------------------------------

# Duplicate setup from Figure 6.7

(nage <- length(growth$age))
# 31
norder.growth <- 6
(nbasis.growth <- nage + norder.growth - 2)
# 35 
(rng.growth <- range(growth$age))
# 1 18 
wbasis.growth <- create.bspline.basis(rng.growth,
                   nbasis.growth, norder.growth,
                   growth$age)

# Specify what to smooth, namely the rate of change of curvature
Lfdobj.growth <- 2 
# Specify smoothing weight (see Figure 6.7 above) 
lambda.gr2.3 <- .03

nage <- length(growth$age)
norder.growth <- 6
nbasis.growth <- nage + norder.growth - 2
rng.growth <- range(growth$age)
# 1 18 
wbasis.growth <- create.bspline.basis(rangeval=rng.growth, norder=norder.growth,
                   breaks=growth$age)

# Smooth all girls;  subset later 
cvec0.growth <- matrix(0,nbasis.growth,1)
Wfd0.growth  <- fd(cvec0.growth, wbasis.growth)
growfdPar2.3 <- fdPar(Wfd0.growth, Lfdobj.growth, lambda.gr2.3)
# Create a functional data object for all the girls
smG.fd.all <- with(growth, smooth.basis(age, hgtf, growfdPar2.3))

# register the girls.  
smGv <- deriv(smG.fd.all$fd, 1)

et.G <- system.time(
smG.reg.1 <- register.fd(smGv,
               WfdParobj=c(Lfdobj=Lfdobj.growth, lambda=lambda.gr2.3))
)
et.G/60

#save(list=c("smG.reg.1", "et.G"), file="registerGirls.Rdata")
#load("registerGirls.Rdata")

class(smG.reg.1)
#  list
sapply(smG.reg.1, class)
#    regfd       Wfd     shift 
#     "fd"      "fd" "numeric"

nPts <- 151
agefine  <- with(growth, seq(age[1], age[nage], len=nPts))

smG.warpmat01 <- eval.monfd(agefine, smG.reg.1$Wfd)

smG.warpmat <- 1+17*smG.warpmat01 / rep(smG.warpmat01[nPts,], each=nPts)
range(smG.warpmat)
#  1 18
str(smG.warpmat)

matplot(agefine, smG.warpmat[, 1:10], type="l")
lines(c(1,18), c(1, 18), lty="dashed", lwd=2, col="black")

# Matches the general look and feel of Figure 6.10
# but not the details.
# This difference is probably due to improvements made
# to the algorithms since those used to produce the book.  

##
## sec. 6.6.  Amplitude and phase variation in growth 
##

# p. 94, Figure 6.11.  Height acceleration for boys and girls
# (cm / year^2)  
smG.regmean <- mean(smG.reg.1$regfd)

# Create a functional data object for all the boys
smB.fd.all <- with(growth, smooth.basis(age, hgtm, growfdPar2.3))

# register the boys
smBv <- deriv(smB.fd.all$fd, 1)

(et.B <- system.time(
smB.reg.1 <- register.fd(smBv,
               WfdParobj=c(Lfdobj=Lfdobj.growth, lambda=lambda.gr2.3))
))
et.B/60

#save(list=c("smB.reg.1", "et.B"), file="registerBoys.Rdata")
#load("registerBoys.Rdata")

smB.regmean <- mean(smB.reg.1$regfd)
str(smB.regmean)
all.equal(smB.regmean$basis, smG.regmean$basis)
# TRUE
all.equal(smB.regmean$fdnames, smG.regmean$fdnames)
# TRUE

smBG.regmean <- with(smB.regmean, fd(cbind(coefs, smG.regmean$coefs),
                                     basis, fdnames) )
plot(smBG.regmean)
plot(deriv(smBG.regmean))

# Not as smooth as Figure 6.11 but similar.
#    A better match might be achieved with more smoothing.
# However, registering the girls took 45 minutes and the boys took 30
# 2007.09.22.  Therefore, I won't push that now.

plot(deriv(smBG.regmean), seq(2, 18, length=101))
# without the initial drop absent from Figure 6.11.  


# p. 95, Figure 6.12.  Registering boys' height velocity to girls

#plot(smBG.regmean)

# After some experimentation, lambda = 0.1,
# gave a rough match to Figure 6.12
smBG.regG.1 <- register.fd(smG.regmean, smBG.regmean, WfdParobj=c(lambda=.1))

sapply(smBG.regG.1, class)
#    regfd       Wfd     shift 
#     "fd"      "fd" "numeric"

smB.warpG.1 <- eval.monfd(agefine, smBG.regG.1$Wfd[1])

smB.wrpG.1 <- 1+17*smB.warpG.1 / rep(smB.warpG.1[nPts,], each=nPts)
range(smB.wrpG.1)
#  1 18
#smG.warpB. <- 1+17*smG.warpB / rep(smG.warpB[nPts,], each=nPts)
#range(smG.warpB.)
#  1 18
str(smB.warpG.1)

op <- par(mfrow=c(2,2))
plot(agefine, smB.wrpG.1, type="l")
lines(c(1,18), c(1, 18), lty="dashed", lwd=2, col="black")

plot(deriv(smBG.regG.1$regfd))

par(op)

# p. 96, Figure 6.13.  Principal components of registered growth acceleration
# for girls

sapply(smG.reg.1, class)

pca.G <- pca.fd(deriv(smG.reg.1$regfd), nharm=3)
class(pca.G)
#[1] "pca.fd"

str(pca.G)
plot(pca.G)

pcaVm.G <- varmx.pca.fd(pca.G)
plot(pcaVm.G)

op <- par(mfrow=c(2,2))
plot(pcaVm.G, cex.main=.9, seq(4, 18, len=101))
par(op)
# Matches Figure 6.13 in broad outline
# but many details are different.
# Might more smoothing produce greater agreement?

# p. 96, Figure 6.14.  PCA of the warping functions

pca.G.w <- pca.fd(deriv(smG.reg.1$Wfd), nharm=3)

pcaVm.G.w <- varmx.pca.fd(pca.G.w)

op <- par(mfrow=c(2,2))
plot(pcaVm.G.w, cex.main=.9, seq(4, 18, len=101))
par(op)

# Differences between the images produced here
# and Figure 6.14 in the book may be due either one of two things:
#
# (a) This script was produced some time after the book was produced,
#     without access to information about exactly what smoothing was used
#     in the book.  More careful experimentation with alternative
#     smoothing might produce a closer match to the book.  
#
# (b) The algorithms have been improved since the book was written,
#     and the current images seem at least as good and perhaps
#     better representations of the data then those in the book.  


