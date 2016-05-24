### R code from vignette source 'deming.Rnw'

###################################################
### code chunk number 1: f1
###################################################
getOption("SweaveHooks")[["fig"]]()
library(deming)
options(continue="  ", width=60)
options(SweaveHooks=list(fig=function() par(mar=c(5.1, 4.1, .3, 1.1))))

tdata <- data.frame(x=1:6, y=c(2.3,1.3, 4.1,3.5,6.3, 3))
lfit <- lm(y ~ x, tdata)
dfit <- deming(y ~ x, tdata)
lfit2 <- lm(x ~ y, tdata)

with(tdata, plot(x, y, xlim=c(0,7), ylim=c(1,7)))
abline(lfit)
abline(-lfit2$coef[1]/lfit2$coef[2], 1/lfit2$coef[2], col=4, lty=2)
abline(dfit, col=2, lwd=2)
segments(tdata$x, tdata$y, tdata$x, predict(lfit), col=1, lty=1)  
segments(tdata$x, tdata$y, predict(lfit2), tdata$y, col=4, lty=2)


###################################################
### code chunk number 2: f2
###################################################
getOption("SweaveHooks")[["fig"]]()
f.ave <- with(ferritin, (old.lot + new.lot)/2)
f.diff<- with(ferritin, old.lot - new.lot)
plot(sqrt(f.ave), f.diff, xaxt='n', 
     xlab="Average", ylab="Difference")
temp <- 0:7*5
axis(1, temp, temp^2)


###################################################
### code chunk number 3: f3
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(f.ave, abs(f.diff/f.ave), log='x', 
     xlab="Average", ylab="Estimated CV")
lines(lowess(f.ave, abs(f.diff/f.ave)), col=2)


###################################################
### code chunk number 4: deming.Rnw:176-186
###################################################
cmat <- matrix(0, nrow=3, ncol=7)
for (i in 1:7) {
    dfit <- deming(new.lot ~ old.lot, data=ferritin, 
                   subset=(period==i), cv=TRUE)
    cmat[1:2,i] <- coef(dfit)
    cmat[3,i] <- coef(lm(new.lot ~ old.lot, ferritin, 
                         subset= (period==i), weight=1/new.lot))[2]
}
dimnames(cmat) <- list(c("Intercept", "old.lot", "old.lot (LS)") , 1:7)
round(cmat,3)


###################################################
### code chunk number 5: deming.Rnw:199-207
###################################################
afit <- deming(aas ~ aes, arsenate, xstd=se.aes, ystd=se.aas)
afit

dfit <- deming(aas ~ aes, arsenate)
lfit <- lm(aas ~ aes, arsenate)
temp <- cbind(coef(afit), coef(dfit), coef(lfit))
dimnames(temp)[[2]] <- c("Ripley", "Deming", "Linear")
round(temp,3)


###################################################
### code chunk number 6: f4
###################################################
getOption("SweaveHooks")[["fig"]]()
tdata <- data.frame(x=1:8, y=c(2,1.2, 4.1,3.5,6.3, 3, 7,6))
xx <- with(tdata, outer(x, x, "-"))
yy <- with(tdata, outer(y, y, "-"))
xx <- c(xx[row(xx) != col(xx)]) #don't compare a point with itself
yy <- c(yy[row(yy) != col(yy)])
plot(xx, yy,
     xlab="Paired x difference", ylab="Paired y difference")
abline(v=0)
fit <- thielsen(y ~ x, tdata)
abline(0, coef(fit)[2], col=2)
text(c(2 ,4, -2, -3), c(5, -4, -5, 4), c("I", "II", "III", "IV"))


###################################################
### code chunk number 7: f5
###################################################
getOption("SweaveHooks")[["fig"]]()
tdata <- data.frame(x=1:8, y=c(2,1.2, 4.1,3.5,6.3, 3, 7,6))
xx <- with(tdata, outer(x, x, "-"))
yy <- with(tdata, outer(y, y, "-"))
xx <- c(xx[row(xx) != col(xx)]) #don't compare a point with itself
yy <- c(yy[row(yy) != col(yy)])
plot(xx, yy,
     xlab="Paired x difference", ylab="Paired y difference")
abline(0, -1)
fit <- pbreg(y ~ x, tdata)
abline(0, coef(fit)[2], col=2)
text(c(0 ,5, 0, -6), c(5, -1, -5, 2), c("I", "II", "III", "IV"))


###################################################
### code chunk number 8: f6
###################################################
getOption("SweaveHooks")[["fig"]]()
zfun <- function(x, y) {
    xx <- outer(x, x, '-')
    yy <- outer(y, y, '-')
    xx <- xx[row(xx) > col(xx)]
    yy <- yy[row(yy) > col(yy)]
    theta <- atan(yy/xx)
    angle <- ifelse(theta >0, theta, theta + pi/2) #in/out angle
    change<- ifelse(theta >0, 1, -1)
    indx <- order(angle)
    list (x= c(0, angle[indx])* 180/pi, 
          y= (sum(theta<0) + c(0, cumsum(change[indx])))/ length(theta))
}
plot(zfun(arsenate$aas, arsenate$aes), type='l',
     xlab="Angle", ylab="Fraction included")
points(zfun(tdata$x, tdata$y))
abline(.5,0)   


###################################################
### code chunk number 9: f7
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(new.lot ~ old.lot, data=ferritin2, subset=(period==2),
     xlab="Old lot", ylab="New lot")
dfit <- deming(new.lot ~ old.lot, ferritin2, subset=(period==2),
               cv=TRUE)

lfit <- lm(new.lot ~ old.lot, ferritin2, subset=(period==2))
pfit <- pbreg(new.lot ~ old.lot, ferritin2, subset=(period==2))
abline(pfit, col=1)
abline(lfit, lty=2)
abline(dfit, lty=3)


