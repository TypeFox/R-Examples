#  ---------------------------------------------------------------------
#            Register the velocity curves for the girls
#  ---------------------------------------------------------------------

#  load the data

load("growthfd")

hgtffdPar <- growthfd$hgtffdPar
hgtffd    <- hgtffdPar$fd

age  <- c( seq(1, 2, 0.25), seq(3, 8, 1), seq(8.5, 18, 0.5))
nage <- length(age)

ncasef <- 54

#  set up the basis for function Wfd

rng     <- c(1,18)
nbasisw <- 15
norder  <- 5
basisw  <- create.bspline.basis(rng, nbasisw, norder)

# set up the mean velocity curve as the preliminary target for
#  registration

hgtfmeanfd <- mean(hgtffd)
y0fd <- deriv(hgtfmeanfd,  1)

#  curves to be registered

yfd  <- deriv(hgtffd, 1)

#  set up functional parameter object for function Wfd

coef0  <- matrix(0,nbasisw,ncasef)
Wfd0   <- fd(coef0, basisw)
lambda <- 10
WfdPar <- fdPar(Wfd0, 2, lambda)

#  register the data.  It might be a good idea to disable
#  buffered output in the Misc menu for the R Console in order
#  to track progress of this fairly slow process.

reglist <- register.fd(y0fd, yfd, WfdPar)

yregfd  <- reglist$regfd  #  registered curves
Wfd     <- reglist$Wfd    #  functions defining warping functions

#  evaluate the registered curves and warping functions

agefine <- seq(1, 18, len=101)
ymat    <- eval.fd(agefine, yfd)
y0vec   <- eval.fd(agefine, y0fd)
yregmat <- eval.fd(agefine, yregfd)
warpmat <- eval.monfd(agefine, Wfd)
warpmat <- 1 + 17*warpmat/(matrix(1,101,1)%*%warpmat[101,])

#  plot the results for each girl:
#    blue:  unregistered curve
#    red:   target curve
#    green: registered curve

par(mfrow=c(1,2),pty="s",ask=T)
for (i in 1:ncasef) {
    plot (agefine, ymat[,i], type="l", ylim=c(0,20), col=4,
          xlab="Year", ylab="Velocity", main=paste("Case",i))
    lines(agefine, y0vec, lty=2, col=2)
    lines(agefine, yregmat[,i],  col=3)
    plot (agefine, warpmat[,i], type="l",
          xlab="Clock year", ylab="Biological Year")
    abline(0,1,lty=2)
}

#  Comments:  we see that not all curves are properly registered.
#     Curves 7, 11, 13 and 25, to mention a few, are so far from
#     the target that the registration is unsuccessful.  This
#     argues for a preliminary landmark registration of the 
#     velocity curves prior to the continuous registration 
#     process.  However, we will see some improvement below.  

#  compute the new mean curve as a target

y0fd2 <- mean(yregfd)

#  plot the unregistered mean and the registered mean

par(mfrow=c(1,1),pty="s",ask=F)
plot(y0fd2, col=4, xlab="Year", ylab="Mean Velocity")
lines(y0fd, col=3)
legend(10,15, c("Registered", "Unregistered"), lty=c(1,1), col=c(4,3))

#  Comment:  The new mean has a sharper peak at the pubertal
#      growth spurt, which is what we wanted to achieve.

#  define the registered curves and the new curves to be registered

yfd2 <- yregfd

#  register the curves again, this time to a better target

reglist2 <- register.fd(y0fd2, yfd2, WfdPar)

yregfd2  <- reglist2$regfd  #  registered curves
Wfd2     <- reglist2$Wfd    #  functions defining warping functions

y0vec2   <- eval.fd(agefine, y0fd2)
yregmat2 <- eval.fd(agefine, yregfd2)
warpmat2 <- eval.monfd(agefine, Wfd2)
warpmat2 <- 1 + 17*warpmat2/(matrix(1,101,1)%*%warpmat2[101,])

#  plot the results for each girl:
#    blue:  unregistered curve
#    red:   target curve
#    green: registered curve

par(mfrow=c(1,2),pty="s",ask=T)
for (i in 1:ncasef) {
    #  plot velocity curves
    plot (agefine, ymat[,i], type="l", ylim=c(0,20), col=4,
          xlab="Year", ylab="Velocity", main=paste("Case",i))
    lines(agefine, y0vec2, lty=2, col=2)
    lines(agefine, yregmat[,i],   col=3, lty=3)
    lines(agefine, yregmat2[,i],  col=3)
    #  plot warping functions
    plot (agefine, warpmat2[,i], type="l",
          xlab="Clock year", ylab="Biological Year")
    abline(0,1,lty=2)
}

#  compute the new mean curve as a target

y0fd3 <- mean(yregfd2)

#  plot the unregistered mean and the registered mean

par(mfrow=c(1,1),pty="s",ask=F)
plot(y0fd3, col=4, xlab="Year", ylab="Mean Velocity")
lines(y0fd2, col=3)
lines(y0fd, col=3, lty=3)
legend(10,15, c("Registered twice", "Registered once", "Unregistered"), 
       lty=c(1,1,3), col=c(4,3,3))

#  Comment:  The second round of registered made hardly any
#    difference for either the individual curves or the mean curve.

