### R code from vignette source 'RobPer_vignette.Rnw'

###################################################
### code chunk number 1: JSS_options_and_package_loading
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("RobPer")


###################################################
### code chunk number 2: Panel_1a
###################################################
data(Mrk421)
temp <- Mrk421
par(mgp=c(2, 1, 0), mar=c(3,3,3,1))
jahre <- seq(1992,2010, by=3)
xlim <-  c(48622, 55000)
yats <- seq(0, 15, by=5)
plot(temp[,1], temp[,2], pch=16, col=rgb(0,0,0,alpha=0.5), axes=F,
    xlab="t (Gregorian date)", ylab="y", xlim=xlim,
    ylim=max(temp[,2]+temp[,3])*c(0, 1.1), cex=0.7)
rect(temp[,1], temp[,2]+temp[,3], temp[,1], temp[,2]-temp[,3],
    border=rgb(0,0,0,alpha=0.5))
mtext("t (Markarian Julian days)", 3, line=2)
axis(2, at=yats)
axis(1, at= as.Date(paste(jahre, "-01-01", sep=""))-
    as.Date("1858-11-17"),labels=jahre)
axis(3)
box()


###################################################
### code chunk number 3: Panel_1b
###################################################
par(mar=c(3,3,0.3,0.1), mgp=c(2,1,0))
p_sam <- 27.31
test <- hist(((Mrk421[,1])%%p_sam), breaks=seq(from=0,
    to=p_sam, length.out=20),main="", col="grey80",freq=F,
    xlim=c(0,p_sam), xlab="Cycle", axes=F)
axis(1, at=(0:9)*3)
axis(2, at=(0:2)*4/100)
y <- test$density
x <- test$mid
fak1 <- sin(x/p_sam*2*pi)
fak2 <- cos(x/p_sam*2*pi)
mo <- lm(y~fak1+fak2)$coeff
ff <- function(x){mo[1]+ mo[2]*sin(x/p_sam*2*pi) + mo[3]*cos(x/p_sam*2*pi)}
curve(ff, add=T, lwd=2)
box()


###################################################
### code chunk number 4: Calculations_for_Figure_2
###################################################
set.seed(923)
zr <- tsgen(ttype="sine", npoints=500, ytype="const", pf=150,
    ncycles=71, ps=4, redpart=0.2, alpha=0,SNR=2, interval=FALSE)
PP <- RobPer(zr, model="sine", regression="L2", weighting=TRUE,
    periods=1:100)
shapes <- betaCvMfit(PP)
myf <- function(x) dbeta(x, shapes[1], shapes[2])
stdf <- function(x) dbeta(x, 1,497/2)


###################################################
### code chunk number 5: Panel_2a
###################################################
par(mar=c(4,4,0.1,0.1))
plot(PP, xlab="Trial period", ylab="Periodogram", type="l", axes=FALSE)
axis(1)
axis(2, at=seq(0,0.06, by=0.02))
abline(h=qbeta(0.95^0.01, 1,497/2), lty=2)
box()


###################################################
### code chunk number 6: Panel_2b
###################################################
par(mar=c(4,4,0.1,0.1))
hist(PP, freq=FALSE, main=" ", xlab="Periodogram",ylab="Density",
    col="grey80", breaks=15, axes=FALSE)
axis(1, at=seq(0,0.06, by=0.02))
axis(2)
box()
curve(myf, add=TRUE)
curve(stdf, add=TRUE, lty=2)


###################################################
### code chunk number 7: Figure_3
###################################################
par(mar=c(3,3,0.3,0.1), mgp=c(2,1,0))
set.seed(12)
PP <- c(rbeta(45, shape1=4, shape2=15), runif(5, min=0.8, max=1))
hist(PP, freq=FALSE, breaks=30, ylim=c(0,7), col="grey90",
    main="", xlab="Periodogram bar")
# true parameters:
myf.true <- function(x) dbeta(x, shape1=4, shape2=15)
curve(myf.true, add=TRUE, lwd=2)
# method of moments:
par.mom <- betaCvMfit(PP, rob=FALSE, CvM=FALSE)
myf.mom <- function(x) dbeta(x, shape1=par.mom[1], shape2=par.mom[2])
curve(myf.mom, add=TRUE, lwd=2, lty=2)
# robust method of moments
par.rob <- betaCvMfit(PP, rob=TRUE, CvM=FALSE)
myf.rob <- function(x) dbeta(x, shape1=par.rob[1], shape2=par.rob[2])
curve(myf.rob, add=TRUE, lwd=2, lty=3)
# CvM distance minimization
par.CvM <- betaCvMfit(PP, rob=TRUE, CvM=TRUE)
myf.CvM <- function(x) dbeta(x, shape1=par.CvM[1], shape2=par.CvM[2])
curve(myf.CvM, add=TRUE, lwd=2, col="grey50")
# Searching for outliers...
abline(v=qbeta((0.95)^(1/50), shape1=par.CvM[1], shape2=par.CvM[2]), col="grey30")
legend("topright", col=c("black", "grey50","black", "black"),lwd=2, lty=c(1,1,3,2),
    legend=c("True", "CvM", "Robust moments", "Moments"))
box()


###################################################
### code chunk number 8: Artificial_example_light_curve
###################################################
library("RobPer")
set.seed(22)
lightcurve <- tsgen(ttype = "sine", ytype = "peak", pf = 7,
    redpart = 0.1, s.outlier.fraction = 0.1, interval = TRUE,
    npoints = 200, ncycles = 25, ps = 20, SNR = 3, alpha = 0)


###################################################
### code chunk number 9: Light_curve_step_by_step
###################################################
set.seed(22)
tt <- sampler(ttype = "sine", npoints = 200, ncycles = 25, ps = 20)


###################################################
### code chunk number 10: periodic_fluctuation
###################################################
yf <- signalgen(tt, ytype = "peak", pf = 7)


###################################################
### code chunk number 11: noise_and_scale_signal
###################################################
temp <- lc_noise(tt, sig = yf, SNR = 3, redpart = 0.1, alpha = 0)
y <- temp$y
s <- temp$s


###################################################
### code chunk number 12: outliers_and_peak
###################################################
temp <- disturber(tt, y, s, ps = 20, s.outlier.fraction = 0.1,
    interval = TRUE)


###################################################
### code chunk number 13: compare_to_tsgen
###################################################
all(cbind(tt, temp$y, temp$s) == lightcurve)


###################################################
### code chunk number 14: Panel_4a
###################################################
par(mar=c(3,3,0.3,0.3), mgp=c(2,1,0))
plot(lightcurve[,1], lightcurve[,2], pch=16, xlab="t", ylab="y",
    main=" ", ylim=range(c(lightcurve[,2]+lightcurve[,3],lightcurve[,2]-lightcurve[,3])),
    col=rgb(0,0,0,alpha=0.5), cex=0.7, axes=FALSE)
axis(1, at=(0:5)*100, labels=c(0, "", 200, "", 400, ""))
axis(2)
rect(lightcurve[,1], lightcurve[,2]+lightcurve[,3], lightcurve[,1],
    lightcurve[,2]-lightcurve[,3], border=rgb(0,0,0,alpha=0.5))
box()


###################################################
### code chunk number 15: Panel_4b
###################################################
par(mar=c(3,3,0.3,0.3), mgp=c(2,1,0))
plot(lightcurve[,1]%%7, lightcurve[,2], pch=16, col=rgb(0,0,0,alpha=0.5),
    xlab="t modulo 7", ylab="y", main=" ",
    ylim=range(c(lightcurve[,2]+lightcurve[,3],lightcurve[,2]-lightcurve[,3])), cex=0.7)
rect(lightcurve[,1]%%7, lightcurve[,2]+lightcurve[,3], lightcurve[,1]%%7,
    lightcurve[,2]-lightcurve[,3], border=rgb(0,0,0,alpha=0.5))


###################################################
### code chunk number 16: Panel_4c
###################################################
par(mar=c(3,3,0.3,0.3), mgp=c(2,1,0))
hist(lightcurve[,1]%%20, xlab="t modulo 20", col="grey", main=" ", freq=FALSE)
dsin <- function(tt) (sin(2*pi*tt/20)+1)/20
curve(dsin, add=TRUE, lwd=2)
box()


###################################################
### code chunk number 17: Periodogram
###################################################
PP <- RobPer(lightcurve, model = "splines", regression = "huber",
    weighting = FALSE, var1 = FALSE, periods = 1:50)


###################################################
### code chunk number 18: Critical_values
###################################################
betavalues <- betaCvMfit(PP)
crit.val <- qbeta((0.95)^(1 / 50), shape1 = betavalues[1],
    shape2 = betavalues[2])


###################################################
### code chunk number 19: Histogram_and_beta_fitting
###################################################
hist(PP, breaks = 20, freq = FALSE, xlim = c(0, 0.08), col = "grey",
    main = "", xlab="Periodogram")
betafun <- function(x) dbeta(x, shape1 = betavalues[1],
    shape2 = betavalues[2])
curve(betafun, add = TRUE, lwd = 2)
abline(v = crit.val, lwd = 2)


###################################################
### code chunk number 20: Method_of_moments_fitting
###################################################
par.mom <- betaCvMfit(PP, rob = FALSE, CvM = FALSE)
myf.mom <- function(x) dbeta(x, shape1 = par.mom[1], shape2 = par.mom[2])
curve(myf.mom, add = TRUE, lwd = 2, lty = 2)
crit.mom <- qbeta((0.95)^(1 / 50), shape1 = par.mom[1],
    shape2 = par.mom[2])
abline(v = crit.mom, lwd = 2, lty = 2)


###################################################
### code chunk number 21: Robust_fitting
###################################################
par.rob <- betaCvMfit(PP, rob = TRUE, CvM = FALSE)
myf.rob <- function(x) dbeta(x, shape1 = par.rob[1], shape2 = par.rob[2])
curve(myf.rob, add = TRUE, lwd = 2, lty = 3)
crit.rob <- qbeta((0.95)^(1 / 50), shape1 = par.rob[1],
    shape2 = par.rob[2])
abline(v = crit.rob, lwd = 2, lty = 3)
legend("topright", lty = 1:3, legend = c("CvM", "Moments",
    "Robust moments"), bg = "white", lwd = 2)
box()


###################################################
### code chunk number 22: Panel_5a
###################################################
par(mar=c(3,3,1,0.7), mgp=c(2,1,0))
hist(PP, breaks = 20, freq = FALSE, xlim = c(0, 0.08), col = "grey",
    main = "", xlab="Periodogram")
curve(betafun, add = TRUE, lwd = 2)
abline(v = crit.val, lwd = 2)
curve(myf.mom, add = TRUE, lwd = 2, lty = 2)
abline(v = crit.mom, lwd = 2, lty = 2)
curve(myf.rob, add = TRUE, lwd = 2, lty = 3)
abline(v = crit.rob, lwd = 2, lty = 3)
legend("topright", lty = 1:3, legend = c("CvM", "Moments", "Robust moments"),
    bg = "white", lwd = 2)
box()


###################################################
### code chunk number 23: Panel_5b
###################################################
par(mar=c(3,3,1,0.1), mgp=c(2,1,0))
plot(1:50, PP, xlab = "Trial period", ylab = "Periodogram", main = "",
    type = "l")
abline(h = crit.val, lwd = 2)
text(7, PP[7]-0.002,7, pos=4)
text(14, PP[14]+0.002,14, pos=4)


###################################################
### code chunk number 24: Periodogram_plot
###################################################
plot(1:50, PP, xlab = "Trial period", ylab = "Periodogram", main = "",
    type = "l")
abline(h = crit.val, lwd = 2)
text(7, PP[7]-0.002,7, pos=4)
text(14, PP[14]+0.002,14, pos=4)


###################################################
### code chunk number 25: Non-robust_periodogram
###################################################
PP <- RobPer(lightcurve, model = "splines", regression = "L2",
    weighting = FALSE, var1 = FALSE, periods = 1:50)


###################################################
### code chunk number 26: Panel_6a
###################################################
betavalues <- betaCvMfit(PP)
crit.val <- qbeta((0.95)^(1 / 50), shape1 = betavalues[1],
    shape2 = betavalues[2])
par(mar=c(3,3,1,0.7), mgp=c(2,1,0))
hist(PP, breaks = 20, freq = FALSE, ylim = c(0, 50),
    col = "grey", main = "", xlab="Periodogram")
betafun <- function(x) dbeta(x, shape1 = betavalues[1],
    shape2 = betavalues[2])
curve(betafun, add = TRUE, lwd = 2)
abline(v = crit.val, lwd = 2)
par.mom <- betaCvMfit(PP, rob = FALSE, CvM = FALSE)
myf.mom <- function(x) dbeta(x, shape1 = par.mom[1],
    shape2 = par.mom[2])
curve(myf.mom, add = TRUE, lwd = 2, lty = 2)
crit.mom <- qbeta((0.95)^(1 / 50), shape1 = par.mom[1],
    shape2 = par.mom[2])
abline(v = crit.mom, lwd = 2, lty = 2)
par.rob <- betaCvMfit(PP, rob = TRUE, CvM = FALSE)
myf.rob <- function(x) dbeta(x, shape1 = par.rob[1],
    shape2 = par.rob[2])
curve(myf.rob, add = TRUE, lwd = 2, lty = 3)
crit.rob <- qbeta((0.95)^(1 / 50), shape1 = par.rob[1],
    shape2 = par.rob[2])
abline(v = crit.rob, lwd = 2, lty = 3)
legend("topright", lty = 1:3, legend = c("CvM", "Moments", "Robust moments"),
    bg = "white", lwd = 2)
box()


###################################################
### code chunk number 27: Panel_6b
###################################################
par(mar=c(3,3,1,0.1), mgp=c(2,1,0))
plot(1:50, PP, xlab = "Trial period", ylab = "Periodogram", main = "",
    type = "l")
abline(h = crit.val, lwd = 2)


###################################################
### code chunk number 28: GROJ0422+32_periodogram
###################################################
data(star_groj0422.32)
PP <- RobPer(star_groj0422.32, periods = 1:330, model = "sine",
    regression = "L2", weighting = FALSE)


###################################################
### code chunk number 29: GROJ0422+32_critical_values
###################################################
shapes <- betaCvMfit(PP)
Crit <- qbeta(0.95^(1 / 330), shape1 = shapes[1], shape2 = shapes[2])


###################################################
### code chunk number 30: NOTE_and_calculations_loading_of_GROJ0422+32
###################################################
#####################
## In order to ease the building of this Vignette, some objects needed
## for the following figures (7 and 8) have previously been calculated
## and can be loaded from grojanalysis.Rdata in the inst/extdata
## directory.
## The calculations can be found in inst/extdata/grojanalysis.R.
#####################
## Objects in grojanalysis.Rdata:
### star_groj0422.32: the light curve
### PP_groj: a list with the three periodograms shown in Figure 7 (b-d)
### Crit_groj: the critical values for Figure 7 (horizontal lines in 7(b-d)
### grojfluc: the light curve with a sine added shown in Figur 8 (a)
### PP_grojfluc: a list with the three periodograms shown in Figure 8 (b-d)
### Crit_grojfluc: the critical values for Figure 8 (horizontal lines in 8(b-d)
## (each (b-d): [[1]] LS-, [[2]] tau-, [[3]] Huber-regression)
#####################
load(system.file("extdata", "grojanalysis.Rdata", package="RobPer"))


###################################################
### code chunk number 31: Panel_7a
###################################################
par(mar=c(3,3,3,1), mgp=c(2,1,0))
xlim_jahr <-  as.numeric(substr(as.Date(range(star_groj0422.32[,1]),
    origin = "1858-11-17", as.character=T),1,4))+c(0,1)
xlim <-  as.Date(paste(xlim_jahr, "-01-01", sep=""))-as.Date("1858-11-17")
plot(star_groj0422.32[,1], star_groj0422.32[,2],  pch=16,
    col=rgb(0,0,0,alpha=0.5),xlab="t (Gregorian date)", ylab="y", xlim=xlim,
    axes=FALSE, cex=0.7)
mtext( "t (Markarian Julian days)", 3, line=2)
axis(1, at= as.Date(paste(xlim_jahr[1]:xlim_jahr[2], "-01-01",  sep=""))-
    as.Date("1858-11-17"),labels=xlim_jahr[1]:xlim_jahr[2])
axis(2)
axis(3)
box()


###################################################
### code chunk number 32: Panel_7b
###################################################
par(mar=c(3,3,1,0.1), mgp=c(2,1,0))
plot(1:330, PP_groj[[1]], ylim=c(0,1.1*Crit_groj[[1]]),
    xlab="Trial period", ylab="Least squares periodogram",
    type="l", cex=0.7)
abline(h=Crit_groj[[1]])


###################################################
### code chunk number 33: Panel_7c
###################################################
par(mar=c(3,3,1,0.1), mgp=c(2,1,0))
plot(1:330, PP_groj[[2]], ylim=c(0,1.1*Crit_groj[[2]]),
    xlab="Trial period", ylab="Tau periodogram", type="l", cex=0.7)
abline(h=Crit_groj[[2]])


###################################################
### code chunk number 34: Panel_7d
###################################################
par(mar=c(3,3,1,0.1), mgp=c(2,1,0))
plot(1:330, PP_groj[[3]],ylim=c(0,1.1*Crit_groj[[3]]),
    xlab="Trial period", ylab="Huber periodogram", type="l", cex=0.7)
abline(h=Crit_groj[[3]])


###################################################
### code chunk number 35: Panel_8a
###################################################
par(mar=c(3,3,3,1), mgp=c(2,1,0))
xlim_jahr <-  as.numeric(substr(as.Date(range(star_groj0422.32[,1]),
    origin = "1858-11-17", as.character=T),1,4))+c(0,1)
xlim <-  as.Date(paste(xlim_jahr, "-01-01", sep=""))-as.Date("1858-11-17")
plot(grojfluc[,1], grojfluc[,2], xlab="t (Gregorian date)", ylab="y",
    pch=16, col=rgb(0,0,0,alpha=0.5), xlim=xlim, axes=FALSE, cex=0.7)
mtext( "t (Markarian Julian days)", 3, line=2)
axis(1, at= as.Date(paste(xlim_jahr[1]:xlim_jahr[2], "-01-01", sep=""))-
    as.Date("1858-11-17"),labels=xlim_jahr[1]:xlim_jahr[2])
axis(2)
axis(3)
box()


###################################################
### code chunk number 36: Panel_8b
###################################################
par(mar=c(3,3,1,0.1), mgp=c(2,1,0))
plot(1:330, PP_grojfluc[[1]], ylim=c(0,Crit_grojfluc[[1]]),
    xlab="Trial period", ylab="Least squares periodogram",
    type="l", cex=0.7)
abline(h=Crit_grojfluc[[1]])


###################################################
### code chunk number 37: Panel_8c
###################################################
par(mar=c(3,3,1,0.1), mgp=c(2,1,0))
plot(1:330, PP_grojfluc[[2]], xlab="Trial period",
    ylab="Tau periodogram", type="l", cex=0.7)
abline(h=Crit_grojfluc[[2]])


###################################################
### code chunk number 38: Panel_8d
###################################################
par(mar=c(3,3,1,0.1), mgp=c(2,1,0))
plot(1:330, PP_grojfluc[[3]], xlab="Trial period",
    ylab="Huber periodogram", type="l", cex=0.7)
abline(h=Crit_grojfluc[[3]])


###################################################
### code chunk number 39: NOTE_and_calculations_loading_of_Mrk421_and_Mrk501
###################################################
#####################
## In order to ease the building of this Vignette, some objects needed
## for the following figures (9 and 10) have previously been calculated
## and can be loaded from MrkAnalysis.Rdata in the inst/extdata
## directory.
## The calculations can be found in inst/extdata/MrkAnalysis.R.
#####################
## Objects in MrkAnalysis.Rdata:
### PPs: list of two lists:
#### 1. with the three periodograms of Mrk421 shown in Figure 9 (a-c)
#### 2. with the three periodograms of Mrk501 shown in Figure 10 (b-d)
### Crits: list of two lists:
#### 1. with the critical values of Mrk421 for Figure 9 (a-c)
#### 2. with the critical values of Mrk501 for Figure 10 (b-d)
### perioden: list of two:
#### 1. trial periods for Mrk421
#### 2. trial periods for Mrk501
## (each (a-c) and (b-d): [[1]] LS-, [[2]] tau-, [[3]] Huber-regression)
#####################
load(system.file("extdata", "MrkAnalysis.Rdata", package="RobPer"))


###################################################
### code chunk number 40: Panel_9a
###################################################
par(mar=c(3,3,0.6,0.5), mgp=c(2,1,0))
plot(perioden[[1]], PPs[[1]][[1]], xlab="Trial period",
    ylab="Least squares periodogram", type="l")
abline(h=Crits[[1]][[1]])


###################################################
### code chunk number 41: Panel_9b
###################################################
par(mar=c(3,3,0.6,0.5), mgp=c(2,1,0))
plot(perioden[[1]], PPs[[1]][[2]], xlab="Trial period",
    ylab="Tau periodogram", type="l")
abline(h=Crits[[1]][[2]])


###################################################
### code chunk number 42: Panel_9c
###################################################
par(mar=c(3,3,0.6,0.5), mgp=c(2,1,0))
plot(perioden[[1]], PPs[[1]][[3]], xlab="Trial period",
    ylab="Huber periodogram", type="l")
abline(h=Crits[[1]][[3]])


###################################################
### code chunk number 43: Mrk501_calculation
###################################################
data(Mrk501)
RobPer(Mrk501, periods = 1:400, model = "step", regression = "L2",
    weighting = FALSE)


###################################################
### code chunk number 44: Panel_10a
###################################################
temp <- Mrk501
par(mgp=c(2, 1, 0), mar=c(3,3,3,0.5))
jahre <- seq(1996,2008, by=1)
xlim <- c(50100,54400)
yats <- seq(0,12, by=3)
plot(temp[,1], temp[,2], pch=16, col=rgb(0,0,0,alpha=0.5), axes=F,
    xlab="t (Gregorian date)", ylab="y", xlim=xlim,
    ylim=max(temp[,2]+temp[,3])*c(0, 1.1), cex=0.7)
rect(temp[,1], temp[,2]+temp[,3], temp[,1], temp[,2]-temp[,3],
    border=rgb(0,0,0,alpha=0.5))
mtext("t (Markarian Julian days)", 3, line=2)
axis(2, at=yats)
axis(1, at= as.Date(paste(jahre, "-01-01", sep=""))-
    as.Date("1858-11-17"),labels=jahre)
axis(3)
box()


###################################################
### code chunk number 45: Panel_10b
###################################################
par(mar=c(3,3,0.6,0.5), mgp=c(2,1,0))
plot(perioden[[2]], PPs[[2]][[1]], xlab="Trial period",
    ylab="Least squares periodogram", type="l", axes=FALSE)
axis(1)
axis(2, at=(0:3)/10)
box()
abline(h=Crits[[2]][[1]])


###################################################
### code chunk number 46: Panel_10c
###################################################
par(mar=c(3,3,0.6,0.5), mgp=c(2,1,0))
plot(perioden[[2]], PPs[[2]][[2]], xlab="Trial period",
    ylab="Tau periodogram", type="l")
abline(h=Crits[[2]][[2]])


###################################################
### code chunk number 47: Panel_10d
###################################################
par(mar=c(3,3,0.6,0.5), mgp=c(2,1,0))
plot(perioden[[2]], PPs[[2]][[3]], xlab="Trial period",
    ylab="Huber periodogram", type="l")
abline(h=Crits[[2]][[3]])


