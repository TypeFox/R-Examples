library(fda)
#  -----------------------------------------------------------------------
#                       Lip Movement Data
#  -----------------------------------------------------------------------
#
#                          Overview of the analyses
#
#  These are rather simple data, involving the movement of the lower lip
#  while saying "bob".  There are 20 replications and 51 sampling points.
#  The data are used to illustrate two techniques:  landmark registration
#  and principal differental analysis.  
#  Principal differential analysis estimates a linear differential equation
#  that can be used to describe not only the observed curves, but also a 
#  certain number of their derivatives.  
#  For a rather more elaborate example of principal differential analysis, 
#  see the handwriting data.
#  -----------------------------------------------------------------------

#  Last modified 2008.06.28;  previously modified 21 March 2006
###
###
### 0.  Access the data:  Instantly available in the 'fda' package 
###
###

###
###
### 1.  Create an 'fd' object 'lipfd'
###
###

##
## 1.1.  Default smooth.basisPar
##
lipfd3 <- smooth.basisPar(liptime, lip)$fd

names(lipfd3$fdnames) <- c("time(seconds)", "replications", "mm")
#op <- par(mfrow=c(2,1), mar=c(5,5,4,2), pty="m", ask=FALSE)
plot(lipfd3,        main="Lip Position", cex=1.2)
plot(lipfd3, Lfd=1, ylab="mm / sec", main="Lip Velocity", cex=1.2)
plot(lipfd3, Lfd=2, ylab="mm / sec / sec", main="Lip Acceleration",
     cex=1.2)
#par(op)

# PROBLEM:  lines too straight, especially position and velocity
# WHY:      Too much smoothing.  
# SOLUTION: Use much less smoothing than the default
##
## 1.2.  Light smoothing
##
lipfd3.12 <- smooth.basisPar(liptime, lip, lambda=1e-12)$fd

names(lipfd3.12$fdnames) <- c("time(seconds)", "replications", "mm")
#op <- par(mfrow=c(2,1), mar=c(5,5,4,2), pty="m", ask=FALSE)
plot(lipfd3.12,        main="Lip Position", cex=1.2)
plot(lipfd3.12, Lfd=1, ylab="mm/sec", main="Lip Velocity", cex=1.2)
plot(lipfd3.12, Lfd=2, ylab="mm/sec/sec", main="Lip Acceleration",
     cex=1.2)
#par(op)

# PROBLEM:  Acceleration not smooth at all ... 
# WHY:      We used cubic splines for location,
#    so the velocity was parabolic splines
#    and acceleration = linear splines (connected straight line segments) 
# SOLUTION: Use quintic splines (degree 5 so order 6) 

##
## 1.3.  Quintic basis (order = 6) 
## 
#lipbasis <- create.bspline.basis(range(liptime), 31, 6) 
#lipfd5 <- smooth.basisPar(liptime, lip, lipbasis, lambda=1e-12)$fd
lipfd5 <- smooth.basisPar(liptime, lip, 6, lambda=1e-12)$fd
names(lipfd5$fdnames) <- c("time(seconds)", "replications", "mm")
#op <- par(mfrow=c(2,1), mar=c(5,5,4,2), pty="m", ask=FALSE)
plot(lipfd5,        main="Lip Position", cex=1.2)
plot(lipfd5, Lfd=1, ylab="mm / sec", main="Lip Velocity", cex=1.2)
plot(lipfd5, Lfd=2, ylab="mm / sec / sec", main="Lip Acceleration",
     cex=1.2)
#par(op)

# PROBLEM:  Acceleration poorly smoothed
# WHY:      The default smoothing operator = int2Lfd(2) = for location
# SOLUTION: Use int2Lfd(4) to smooth acceleration of acceleration  

##
## 1.4.  Penalize the 4th derivative, not the second
##
lipfd <- smooth.basisPar(liptime, lip, 6, Lfdobj=int2Lfd(4),
                         lambda=1e-12)$fd
names(lipfd$fdnames) <- c("time(seconds)", "replications", "mm")
#op <- par(mfrow=c(2,1), mar=c(5,5,4,2), pty="m", ask=FALSE)
plot(lipfd,        main="Lip Position", cex=1.2)
plot(lipfd, Lfd=1, ylab="mm / sec", main="Lip Velocity", cex=1.2)
plot(lipfd, Lfd=2, ylab="mm / sec / sec", main="Lip Acceleration",
     cex=1.2)
#par(op)

##
## 1.5.  plotfit.fd?
##
plotfit.fd(lip, liptime, lipfd)

plotfit.fd(lip, liptime, lipfd, residual=TRUE, type='b',
           sortwrd=TRUE)

##
## 2.  Register the data
##
#  --------------------------------------------------------------------
#       Register the data using the two landmarks defined by
#        the left and right elbows.  
#  --------------------------------------------------------------------

# Optionally:  Manually identify these points in each curve

#par(mfrow=c(1,1),pty="m")
#lipmarks <- matrix(0,20,nmarks)
#index <- 1:20
#for (i in index) {
#  plot(liptime, lipmat[,i], xlab="", ylab="", main=paste("Curve",i))
#  indexi <- identify(liptime, lipmat[,i], n=nmarks)
#  lipmarks[i,] <- liptime[indexi]
#}

lipmeanmarks <- apply(lipmarks,2,mean)

#  -------------   register the curves  --------------------

#  First create a basis object for the warping function
#  it has order 4 (piecewise cubic) and two interior knots
#  positioned at the mean landmark values since
#  NBASIS = NORDER + # interior knots

wnbasis <- 6
wnorder <- 4
wbreaks <- c(0,lipmeanmarks,0.35) 

#warpbasis <- create.bspline.basis(liprange, wnbasis, wnorder, wbreaks);
#warpbasis <- create.bspline.basis(range(lip), wnbasis, wnorder, wbreaks);
warpbasis <- create.bspline.basis(nbasis=wnbasis, norder=wnorder,
                                  breaks=wbreaks);
fd(basisobj=warpbasis)
WfdPar    <- fdPar(fd(basisobj=warpbasis), 2, 1e-4)
WfdPar.    <- fdPar(fd(matrix(0,wnbasis,1), warpbasis), 2, 1e-4)
all.equal(WfdPar, WfdPar.)

lipreglist <- landmarkreg(lipfd, as.matrix(lipmarks), lipmeanmarks, WfdPar)

lipregfd   <- lipreglist$regfd
lipwarpfd  <- lipreglist$warpfd

#  plot unregistered and registered curves

par(mfrow=c(1,2), pty="s")

plot(lipfd, main="Unregistered")
lines.fd(lipmeanfd, lty=2)
abline(v=lipmeanmarks,lty=2)

plot(lipregfd, main="Registered")
lines.fd(lipmeanfd, lty=2)
abline(v=lipmeanmarks,lty=2)

#  plot warping functions and deformations

par(mfrow=c(1,2), pty="s")
plot(lipwarpfd, href=FALSE, main="Warping Functions")
abline(0,1,lty=2)
hmat <- eval.fd(liptime, lipwarpfd)
defmat <- hmat - outer(liptime,rep(1,20))
matplot(liptime,defmat,type="l",lty=1, 
        xlab="Normalized time", ylab="Warped Normalized time",
        main="Deformation Functions")
abline(h=0,lty=2)
##
## 3.  Principal Components Analysis
##
#  ------------  carry out a pca and plot results  -------------------

lambda    <- 1e-6
pcafdPar  <- fdPar(lipbasis, 2, lambda)
lippca.fd <- pca.fd(lipfd, nharm=3, pcafdPar)

par(mfrow=c(1,1),pty="m")
plot.pca.fd(lippca.fd)

lipeigvals <- lippca.fd[[2]]
plot(1:19, log10(lipeigvals[1:19]), type="b",
     xlab="Eigenvalue Number", ylab="", main="Log10 Eigenvalues")
##
## 4.  Principal Differential Analysis
## 
#  ---------------------------------------------------------------------
#                    Principal differential analysis  
#  ---------------------------------------------------------------------

#  set up a second order linear differnetial equation solution

liprange = range(liptime)

pdabasisfd <- create.bspline.basis(liprange, nbasis=21)
betafdPar  <- fdPar(pdabasisfd)

#  set up list of functional parameter objects for weight fns.

bwtlist = vector("list", 2)
bwtlist[[1]] <- betafdPar
bwtlist[[2]] <- betafdPar

xfdlist <- list(lipfd)

pdaList <- pda.fd(xfdlist, bwtlist)

#  plot weight functions

bwtestlist <- pdaList$bwtlist

par(mfrow=c(2,1),pty="m")
for (j in 1:2) {
	bfdParj <- bwtestlist[[j]]
	bvals = eval.fd(liptime,bwtestlist[[j]]$fd)
	plot(liptime,bvals,type='l')
}

#  compute forcing functions

Lfdest <- Lfd(2, bwtestlist)

force        <- eval.fd(liptime, lipfd, Lfdest)
lipaccel     <- eval.fd(liptime, lipfd, 2)
lipmeanaccel <- apply(lipaccel, 1, mean)

par(mfrow=c(1,1),ask=FALSE)
yrange <- c(min(min(lipmeanaccel),min(force)),
            max(max(lipmeanaccel),max(force)))
matplot(liptime, force, type="l", lty=1, ylim=yrange)
lines(liptime, lipmeanaccel, lty=4, lwd=2)

#  plot the mean forcing function along with second deriv.

forcemean <- apply(force, 1, mean)

plot(liptime, forcemean, type="l", lty=1, ylim=yrange)
lines(liptime, lipmeanaccel, lty=4)

#  solve equation

result <- odesolv(bwtestlist)
xp <- result[[1]]
yp <- result[[2]]

#  plot the two solutions

par(mfrow=c(2,1),pty="m")
pltrng <- c(min(yp[1,,]), max(yp[1,,]))
matplot(xp,t(yp[1,,]), type="l", lty=1, ylim=pltrng, main="Function")
abline(h=0, lty=2)
pltrng <- c(min(yp[2,,]), max(yp[2,,]))
matplot(xp,t(yp[2,,]), type="l", lty=1, ylim=pltrng, main="Derivative")
abline(h=0, lty=2)

#  plot fit to each curve

lipmat   <- eval.fd(liptime, lipfd)
D2lipmat <- eval.fd(liptime, lipfd, 2)

umat <- matrix(0,length(liptime),2)
umat[,1] <- approx(xp, t(yp[1,1,]), liptime)$y
umat[,2] <- approx(xp, t(yp[1,2,]), liptime)$y

par(mfrow=c(1,2),pty="s",ask=TRUE)
index <- 1:20
for (i in index) {
    plot(liptime, force[,i], type="l",
          ylim=c(-1000,1000), xlab="Normalized Time", ylab="", 
          main=paste("Record",i,"Forcing Fn."))
    lines(liptime, D2lipmat[,i],lty=4)
    abline(h=0,lty=2)
    xhat <- lipmat[,i] - lsfit(umat, lipmat[,i], int=FALSE)$residual
    matplot(liptime, cbind(xhat, lipmat[,i]), type="l", lty=c(1,2),
          xlab="Normalized Time", ylab="", main="Function")
}

