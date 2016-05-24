### R code from vignette source 'seewave_analysis.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: seewave_analysis.rnw:77-79
###################################################
options(SweaveHooks=list(fig=function()
par(mar=c(5.1, 4.1, 1, 1))))


###################################################
### code chunk number 2: samplingfrequency
###################################################
getOption("SweaveHooks")[["fig"]]()
library(seewave)
s.44kHz <- synth(d=0.005, cf=440, f=44100, out="Wave")
s.11kHz <- synth(d=0.005, cf=440, f=22050, out="Wave")
layout(matrix(1:2, nr=2), height=c(0.8,1))
par(mar=c(0, 4, 1, 1), oma=c(1,1,0,0))
oscillo(s.44kHz, type="o", xaxt="n", tlab="", alab="", cex=0.5)
par(mar=c(4.5,4,0,1))
oscillo(s.11kHz, type="o", tlab="", alab="", cex=0.5)
mtext("Time (s)", side=1, out=TRUE, at=0.55, line=-1)
mtext("Amplitude", side=2, out=TRUE, line=-2, at=0.5)


###################################################
### code chunk number 3: quantisation
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar=c(1,4,1,1))
plot(sin(seq(0,2*pi,length=100)), type="l", ann=FALSE, xaxt="n", yaxt="n", col="red", lwd=2, ylim=c(-1.3,1))
par(new=TRUE)
barplot(c(c(1,2,3,3,3,3,2,1)+0.5,c(-1,-1,-2,-3,-3,-3,-3,-2,-1)-.5),yaxt="n", space=1.5, ylim=c(-4.5,3.5), col="#00000025")
y <- seq(-4.5, 3.5, by=1)
#y <- seq(-3,3, len=2^3+1)
abline(h=y, col="grey")
abline(h=0, lwd=2)
axis(side=2, at=y[-length(y)]+0.5, label=c("100", "101", "110", "111", "000", "001", "010", "011"), las=1, tick=FALSE)
box()


###################################################
### code chunk number 4: envelope
###################################################
getOption("SweaveHooks")[["fig"]]()
data(tico) 
env(tico, envt="abs", colwave="grey")
par(new=TRUE)
env(tico, colwave=2)
legend("topright", lty=1, col=c("grey",2), legend=c("Absolute", "Hilbert"), bty="n")


###################################################
### code chunk number 5: timer
###################################################
getOption("SweaveHooks")[["fig"]]()
timer(tico,f=22050,threshold=5,msmooth=c(50,0))


###################################################
### code chunk number 6: s
###################################################
getOption("SweaveHooks")[["fig"]]()
d <- 0.05
f <- 8000
s1 <- synth(d=d, f=f, cf=1000, a=1, out="Wave")
s2 <- synth(d=d, f=f, cf=2000, a=0.5, out="Wave")
s3 <- synth(d=d, f=f, cf=3000, a=0.25, out="Wave")
s <- s1+s2+s3
oscillo(s)


###################################################
### code chunk number 7: s-decomposed
###################################################
getOption("SweaveHooks")[["fig"]]()
op <- layout(matrix(1:3, nr=3), height=c(1,0.8,1.2))
par(mar=c(0,4,1,1), oma=c(1,0,0,0))
oscillo(s1, labels=FALSE,xaxt="n", bty="o")
legend("topright", legend="1 kHz", bty="n", bg="white", cex=3, text.col=2)
par(mar=c(0,4,0,1))
oscillo(s2, labels=FALSE,xaxt="n", bty="o")
legend("topright", legend="2 kHz", bty="n", bg="white", cex=3, text.col=2)
par(mar=c(4.5,4,0,1))
oscillo(s3, labels=FALSE, bty="o")
legend("topright", legend="3 kHz", bty="n", bg="white", cex=3, text.col=2)
mtext("Time (s)", side=1, out=TRUE, at=0.55, line=-1)
mtext("Amplitude", side=2, out=TRUE, line=-2, at=0.5)
par(op)


###################################################
### code chunk number 8: s-spec
###################################################
getOption("SweaveHooks")[["fig"]]()
x <- 1:3
y <- c(1, 0.5, 0.25)
plot(x, y, xlim=c(0,4), ylim=c(0,1), pch=19, xlab=expression(paste("Frequency ", omega[i], " (kHz)")), ylab=expression(paste("Relative amplitude ", a[i], " (no scale)")), las=1)
segments(x0=x,y0=0, y1=y)


###################################################
### code chunk number 9: DTFT
###################################################
getOption("SweaveHooks")[["fig"]]()
spec(tico, col="grey")
local <- spec(tico, wl = 512, at = 1.1, plot=FALSE)
mean <- meanspec(tico, plot=FALSE)
lines(local, col=2)
lines(mean, col=4)
legend("topright", legend=c("DTFT on complete sound", "DTFT on a sound section", "Mean spectrum (STFT)"), lty=1, col=c("grey", "red", "blue"), bty="n")


###################################################
### code chunk number 10: ftwindow
###################################################
getOption("SweaveHooks")[["fig"]]()
a<-ftwindow(512,wn="rectangle")
b<-ftwindow(512,wn="bartlett")
c<-ftwindow(512, wn="hanning")
all<-cbind(a,b,c)
matplot(all,type="l", col=c("green", "red", "blue"), lty=1, ylab="Amplitude", xlab="Index")
legend(x=380,y=0.95, legend=c("Rectangle", "Bartlett", "Hamming"), col=c("green", "red", "blue"), lty=1, bty="n")


###################################################
### code chunk number 11: ftwindow-test
###################################################
getOption("SweaveHooks")[["fig"]]()
dB <- "max0"
spec(s, wn="rectangle", dB=dB, alim=c(-100,0), col="darkgreen")
bartlett <- spec(s, wn="bartlett", dB=dB, plot=FALSE)
lines(bartlett, col="red")
hamming <- spec(s, wn="hanning", dB=dB, plot=FALSE)
lines(hamming, col="blue")
legend("topright", legend=c("Rectangle", "Bartlett", "Hanning"), col=c("darkgreen", "red", "blue"), lty=1, bty="n", cex=0.75)


###################################################
### code chunk number 12: dynspec (eval = FALSE)
###################################################
## dynspec(tico, wl=1024, osc=TRUE)


###################################################
### code chunk number 13: spectro-grid
###################################################
getOption("SweaveHooks")[["fig"]]()
spectro <- spectro(tico, grid=FALSE, scale=FALSE)
abline(v=seq(0, 1.795, length=78), col="lightgrey")
box()


###################################################
### code chunk number 14: spectro-grid-ovlp
###################################################
getOption("SweaveHooks")[["fig"]]()
spectro.ovlp <- spectro(tico, ovlp=50, scale=FALSE, grid=FALSE)
abline(v=seq(0, 1.795, length=78*2), col="lightgrey")
box()


###################################################
### code chunk number 15: spectro-components
###################################################
getOption("SweaveHooks")[["fig"]]()
op <- par(mar=rep(1,4))
filled.contour(x=spectro$time, y=spectro$freq, z=t(spectro$amp))
par(op)


###################################################
### code chunk number 16: spectro-spec
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(x=spectro$freq, y=spectro$amp[,15], type="o", xlab="Frequency (kHz)", ylab="Relative amplitude (dB)")


###################################################
### code chunk number 17: zc
###################################################
getOption("SweaveHooks")[["fig"]]()
s <- synth(d=0.05, f=8000, cf=440)
s.interpol <- approx(s, n = length(s)*10)
s.interpol <- as.matrix(s.interpol$y)
layout(matrix(1:2, nr=2), height=c(0.8,1))
par(mar=c(0,4,1,1))
oscillo(s, f = 8000, type="o", cex = 0.75, to=0.01, xaxt="n", alab="", tlab="")
arrows(x0=0.0046, x1=0.006855, y0=0, y1=0, col="red", angle=90, length=0.1, lwd=1.5, code=3)
par(mar=c(4.5,4,0,1))
oscillo(s.interpol, f = 8000*10,  cex = 0.75,  type="o", to=0.01, alab="", tlab="")
arrows(x0=0.0046, x1=0.006855, y0=0, y1=0, col="red", angle=90, length=0.1, lwd=1.5, code=3)
mtext("Time (s)", side=1, out=TRUE, at=0.55, line=-1)
mtext("Amplitude", side=2, out=TRUE, line=-2, at=0.6)


###################################################
### code chunk number 18: zc-tico
###################################################
getOption("SweaveHooks")[["fig"]]()
f <- 22050
sel <- cutw(tico, f=f, from=0.5, to=0.9)
layout(matrix(1:2, nr=2), height=c(0.8,1))
par(mar=c(0, 4, 1, 1), oma=c(1,1,0,0))
zc(sel, f=f, threshold=5, xlab="", ylab="", xaxt="n", cex=0.75, warning=FALSE)
par(mar=c(4.5, 4, 0, 1))
zc(sel, f=f, threshold=5, interpol=20,xlab="", ylab="", cex=0.75, warning=FALSE)
mtext("Time (s)", side=1, out=TRUE, at=0.55, line=-1)
mtext("Frequency (kHz)", side=2, out=TRUE, line=-1, at=0.6)


###################################################
### code chunk number 19: ceps
###################################################
getOption("SweaveHooks")[["fig"]]()
data(peewit)
par(mar=c(4.5,4.5,4.5,2))
ceps(peewit, at=0.4, wl=1024, col=2, alim=c(0,2000))


