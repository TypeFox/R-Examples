### R code from vignette source 'DDHFm.Rnw'

###################################################
### code chunk number 1: lib
###################################################
library("DDHFm")
data(cdna)


###################################################
### code chunk number 2: cdna-raw
###################################################
cdna.mn <- apply(cdna,1,mean)
cdna.sd <- apply(cdna,1,sd)
plot(cdna.mn, cdna.sd, xlab="Replicates mean", ylab="Replicates sd")


###################################################
### code chunk number 3: cdna-ddhfm (eval = FALSE)
###################################################
## #
## # Do the DDHFm transform
## #
## cdna.dd <- DDHFm(cdna)
## #
## # Set some graphics parameters 
## #
## op <- par(fg="gray90", tcl=-0.2, mgp=c(3,0.5,0))
## #
## # Compute replicate means and standard deviations
## #
## cdna.dd.mn <- apply(cdna.dd,1,mean)
## cdna.dd.sd <- apply(cdna.dd,1,sd)
## #
## # Plot the rep means and sds
## #
## plot(cdna.dd.mn, cdna.dd.sd, xlab="DDHFm rep mean", ylab="DDHFm rep sd",
## 	col="gray90")
## library("lokern")
## cdna.dd.lo <- lokerns(x=cdna.dd.mn, y=cdna.dd.sd)
## lines(cdna.dd.lo$x.out, cdna.dd.lo$est, col=2)


###################################################
### code chunk number 4: cdna-vsn (eval = FALSE)
###################################################
## library("vsn")
## cdna.vsn <- vsn(cdna)
## #
## # Compute replicate means and standard deviations
## #
## cdna.vsn.mn <- apply(cdna.vsn,1,mean)
## cdna.vsn.sd <- apply(cdna.vsn,1,sd)
## #
## # Plot the rep means and sds
## #
## plot(cdna.vsn.mn, cdna.vsn.sd, xlab="vsn rep mean", ylab="vsn rep sd",
## 	col="gray90")
## cdna.vsn.lo <- lokerns(x=cdna.vsn.mn, y=cdna.vsn.sd)
## lines(cdna.vsn.lo$x.out, cdna.vsn.lo$est, col=3)


###################################################
### code chunk number 5: cdna-vsn-compare (eval = FALSE)
###################################################
## cdna.vsn.lo.withsamebandwidthasdd <- lokerns(x=cdna.vsn.mn, y=cdna.vsn.sd,
## 		inputb=TRUE, bandwidth=cdna.dd.lo$bandwidth)
## plot(seq(from=0, to=1, length=300), cdna.vsn.lo$est, ylim=c(0,1.1),
## 	xlab="x", ylab="VSN transformed scale", type="n")
## lines(seq(from=0, to=1, length=300), cdna.vsn.lo$est, col=3)
## lines(seq(from=0, to=1, length=300),
## 	cdna.vsn.lo.withsamebandwidthasdd$est, col=4)
## lines(seq(from=0, to=1, length=300), cdna.dd.lo$est*0.4988, col=2)


###################################################
### code chunk number 6: pivmkfn
###################################################
makepiv <- function(nbit=4){
p1 <- rep(3, nbit)
p2 <- rep(5,nbit)
p3 <- rep(1,nbit)
p4 <- rep(10, nbit)
xx <- seq(from=1, to=10, length=nbit*4)/3
p5  <- xx^2
c(p1,p2,p3,p4,p5)
}


###################################################
### code chunk number 7: pivmake
###################################################
piv <- makepiv(nbit = 32)
l <- length(piv)
ts.plot(piv, xlab = "x", ylab = "Intensity vector")


###################################################
### code chunk number 8: pseq
###################################################
pseq <- rpois(l, lambda = piv)
plot(1:l, pseq, xlab = "x", ylab = "Poisson sample from intensity vector")
lines(1:l, piv, lty = 2)


###################################################
### code chunk number 9: pddhf
###################################################
pseq.ddhf <- ddhft.np.2(pseq)
plot(pseq.ddhf$mu, pseq.ddhf$sigma2, xlab="Estimated mu", ylab="Estimated sd")
lines(pseq.ddhf$mu, pseq.ddhf$sigma)
lines(pseq.ddhf$mu, sqrt(pseq.ddhf$mu), lty = 2)


###################################################
### code chunk number 10: psmooth
###################################################
pseq2.ddhf <- pseq.ddhf
library("wavethresh")
hftwd <- wd(pseq.ddhf$hft, filter.number = 1, family = "DaubExPhase")
madmad <- function(x) mad(x)^2
hftwdT <- threshold(hftwd, policy = "universal", levels = hftwd$nlevels - 1, dev = madmad, return.thresh = TRUE)
hftwd.thresh <- threshold(hftwd, policy = "manual", value = hftwdT)
hftwr <- wr(hftwd.thresh)
pseq2.ddhf$hft <- hftwr
plot(1:l, pseq.ddhf$hft, xlab="x", ylab="Haar wavelet fit to DDHF data")
lines(1:l, hftwr)


###################################################
### code chunk number 11: ptrudomain
###################################################
pest2 <- ddhft.np.inv(pseq2.ddhf)
plot(1:l, pseq, xlab = "x", ylab = "Poisson data")
lines(1:l, pest2)
lines(1:l, piv, col = 2, lty = 2)
lines(1:l, pss <- smooth.spline(1:l, y = pseq)$y, col = 3)
hfssq <- sum((pest2 - piv)^2)
ssssq <- sum((pss - piv)^2)
title(sub = paste("SSQ: DDHF=", round(hfssq, 0), " SS=", 
        round(ssssq, 0)))


