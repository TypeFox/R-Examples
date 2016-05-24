### R code from vignette source 'mod2user.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: mod2user.Rnw:34-38
###################################################
require(lmodel2)
options(width=72)
figset <- function() par(mar=c(4,4,1,1)+.1)
options(SweaveHooks = list(fig = figset))


###################################################
### code chunk number 2: mod2user.Rnw:470-473
###################################################
data(mod2ex1)
Ex1.res <- lmodel2(Predicted_by_model ~ Survival, data=mod2ex1, nperm=99)
Ex1.res


###################################################
### code chunk number 3: mod2user.Rnw:476-479
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(Ex1.res, centro = TRUE, xlab="log(observed survival time", ylab="Forecasted")
abline(diff(colMeans(mod2ex1)), 1, lty=2)
legend("topleft", c("MA regression", "Confidence limits", "45 degree line"), col=c("red", "grey", "black"), lty=c(1,1,2))


###################################################
### code chunk number 4: mod2user.Rnw:535-538
###################################################
data(mod2ex2)
Ex2.res = lmodel2(Prey ~ Predators, data=mod2ex2, "relative","relative",99)
Ex2.res


###################################################
### code chunk number 5: mod2user.Rnw:541-546
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(Ex2.res, confid=FALSE, xlab="Eagle rays (predators)", ylab="Bivalves (prey)", main = "", centr=TRUE)
lines(Ex2.res, "OLS", col=1, confid=FALSE)
lines(Ex2.res, "SMA", col=3, confid=FALSE)
lines(Ex2.res, "RMA", col=4, confid=FALSE)
legend("topleft", c("OLS", "MA", "SMA", "RMA"), col=1:4, lty=1)


###################################################
### code chunk number 6: mod2user.Rnw:622-625
###################################################
data(mod2ex3)
Ex3.res = lmodel2(No_eggs ~ Mass, data=mod2ex3, "relative", "relative", 99)
Ex3.res


###################################################
### code chunk number 7: mod2user.Rnw:628-632
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(Ex3.res, method="OLS", conf = FALSE, centroid=TRUE, main="", xlab = "Fish mass (x100 g)", ylab = "No. of eggs", col=1)
lines(Ex3.res, "SMA", col=2, conf = FALSE)
lines(Ex3.res, "RMA", col=3, conf = FALSE)
legend("topleft", c("OLS", "SMA", "RMA"), col=1:3, lty=1)


###################################################
### code chunk number 8: mod2user.Rnw:664-667
###################################################
data(mod2ex4)
Ex4.res = lmodel2(y ~ x, data=mod2ex4, "interval", "interval", 99)
Ex4.res


###################################################
### code chunk number 9: mod2user.Rnw:703-706
###################################################
data(mod2ex5)
Ex5.res <- lmodel2(random_y ~ random_x, data=mod2ex5,"interval","interval",99)
Ex5.res


###################################################
### code chunk number 10: mod2user.Rnw:709-714
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(Ex5.res, conf = FALSE, main = "", centro = TRUE)
lines(Ex5.res, "OLS", conf=F, col=1)
lines(Ex5.res, "RMA", conf=F, col=4)
lines(Ex5.res, "SMA", conf=F, col=3)
legend("topright", c("OLS","MA","SMA", "RMA"), col=1:4, lty=1)


