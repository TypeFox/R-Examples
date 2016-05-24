# ---------------------------------------------------------------
#                   Figure Refinery
# ---------------------------------------------------------------

# Last modified 2008.12.26 by Spencer Graves
# Previously modified  5 November 2008

#  attach the FDA functions

# 'refinery' is a data set in 'fda' so no need to read it.
#refinery <- t(matrix(scan("../data/refinery.txt", 0), 3, 194))

tval <- refinery[,1]  #  observation time
uval <- refinery[,2]  #  reflux flow
yval <- refinery[,3]  #  tray 47 level
n <- length(tval)

#  center the data on mean values prior to change

uval <- uval - mean(uval[1:60])
yval <- yval - mean(yval[1:60])

#  plot data

op <- par(mfrow=c(2,1), pty="m")
plot(tval, yval, type="p", ylim=c(-0.5,4.5),
     xlab="Time", ylab="Tray 47 level")

plot(tval, uval, type="p", ylim=c(-0.6,0.2),
     xlab="Time", ylab="Reflux flow")
par(op)

#  ------------------------------------------------------------
#    A forced first order constant coefficient DIFE
#        exp(-wt)[x.0 + D^{-1} [exp(wt) u(t)]
#    Annihilating operator is L = w x + Dx
#    The data are the Corpus Christi refinery data
#  ------------------------------------------------------------

refrange <- c(tval[1], tval[n])
delta    <- 1/(n-1)
tbreak   <- tval[67]

#  set up basis for input variable

norder  <- 1
nubasis <- 2
ubreaks <- c(refrange[1], tbreak, refrange[2])
ubasis  <- create.bspline.basis(refrange, nubasis, norder, ubreaks)

ufd <- smooth.basis(tval, uval, ubasis)$fd

#  set up basis for the output variable
#  put three coincident knots at tbreak

norder  <-  4
yknots  <- c(seq(0, tbreak, len=3), tbreak,
             seq(tbreak, refrange[2], len=5))
nybasis <- length(yknots) + norder - 2
ybasis  <- create.bspline.basis(refrange, nybasis, norder, yknots)

yfd <- smooth.basis(tval, yval, ybasis)$fd

yvec <- eval.fd(tval, yfd)

#  plot the data with fits and knots

par(mfrow=c(1,1))
plot(tval, yval, type="p", xlab="minutes",
     ylab="Tray 47 level", xlim=refrange)
lines(tval, yvec, lwd=3)
for (i in 2:length(yknots)-1)
    lines(c(yknots[i],yknots[i]), c(0,4), lty=4, lwd=2)


#  set up fRegress concurrent model analysis,
#  with the dependent variable being tray 47 level
#  and independent variable being reflux flow

p <- 1

norder    <- 4
bbreaks   <- c(seq(0, tbreak, len=2), tbreak,
               seq(tbreak, refrange[2], len=5))
nbbasis   <- length(bbreaks) + norder - 2
bbasis    <- create.bspline.basis(refrange, nbbasis, norder, bbreaks)
betafd    <- fd(matrix(0,nbbasis,p), bbasis)
betafdPar <- fdPar(betafd, 2, 1e-4)

betalist <- list(betafdPar)

xfdlist  <- list(ufd)

#  carry out the analysis

fRegressList <- fRegress(yfd, xfdlist, betalist)

betaestlist <- fRegressList$betaestlist
yhatfdobj   <- fRegressList$yhatfdobj

#  plot the regression function

betafdPar <- betaestlist[[1]]
betafd    <- betafdPar$fd

plot(betafd, xlab="minutes", ylab="beta(t)")

#  set up fRegress concurrent model analysis,
#  with the dependent variable being derivative of tray 47 level
#  and independent variables being tray 47 level and reflux flow

nDybasis <- length(yknots) + 3 - 2
Dybasis  <- create.bspline.basis(refrange, nDybasis, 3, yknots)

tvalD <- c(tval[1:66], tbreak - 1e-5, tbreak + 1e-5, tval[68:n])
Dyvec <- eval.fd(tvalD, yfd, 1)
Dyfd  <- smooth.basis(tvalD, Dyvec, Dybasis)$fd
Dyfdhat <- eval.fd(tvalD, Dyfd)

p <- 2

bbasis    <- create.constant.basis(refrange)
betafd    <- fd(0, bbasis)
betafdPar <- fdPar(betafd)

betalist[[1]] <- betafdPar
betalist[[2]] <- betafdPar

xfdlist[[1]] <- yfd
xfdlist[[2]] <- ufd

#  carry out the analysis

fRegressList <- fRegress(Dyfd, xfdlist, betalist)

betaestlist <- fRegressList$betaestlist
Dyhatfdobj  <- fRegressList$yhatfdobj

betafdPar1 <- betaestlist[[1]]
betafdPar2 <- betaestlist[[2]]
beta1 <- betafdPar1$fd$coefs
beta2 <- betafdPar2$fd$coefs

print(paste("Beta_1 =",round(beta1,3)))
print(paste("Beta_2 =",round(beta2,3)))

#  plot the first derivative and fit

#Dyhatvec <- eval.fd(tvalD, Dyhatfdobj)
Dyhatvec <- predict(Dyhatfdobj, tvalD)

plot(Dyfd, ylim=c(-.02,0.1))
lines(tvalD, Dyhatvec, lwd=3)

#  plot the data and the fit

nD <- length(tvalD)
yhatvec <- rep(0,nD)
yhatvec[68:nD] <- 0.4924*(beta2/beta1)*(1-exp(beta1*(tval[68:nD]-67)))

#  plot the data and fit

par(mfrow=c(1,1))
plot(tval, yval, type="p", xlim=refrange,
     xlab="minutes", ylab="Tray 47 level")
lines(tvalD, yhatvec, lwd=3)


