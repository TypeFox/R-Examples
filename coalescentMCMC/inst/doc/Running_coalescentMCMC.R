### R code from vignette source 'Running_coalescentMCMC.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Running_coalescentMCMC.Rnw:19-20
###################################################
options(width=60)


###################################################
### code chunk number 2: Running_coalescentMCMC.Rnw:57-59
###################################################
library(coalescentMCMC)
args(coalescentMCMC)


###################################################
### code chunk number 3: Running_coalescentMCMC.Rnw:81-83
###################################################
data(woodmouse)
out <- coalescentMCMC(woodmouse, ntrees = 300, burnin = 100, printevery = 0)


###################################################
### code chunk number 4: Running_coalescentMCMC.Rnw:106-107
###################################################
plot(out)


###################################################
### code chunk number 5: Running_coalescentMCMC.Rnw:112-114
###################################################
TR <- getMCMCtrees()
TR


###################################################
### code chunk number 6: Running_coalescentMCMC.Rnw:118-120
###################################################
dim(out)
colnames(out)


###################################################
### code chunk number 7: Running_coalescentMCMC.Rnw:125-126
###################################################
out2 <- coalescentMCMC(woodmouse, ntrees = 300, burnin = 100, model = "time", printevery = 0)


###################################################
### code chunk number 8: Running_coalescentMCMC.Rnw:145-146
###################################################
plot(out2)


###################################################
### code chunk number 9: Running_coalescentMCMC.Rnw:149-151
###################################################
dim(out2)
colnames(out2)


###################################################
### code chunk number 10: Running_coalescentMCMC.Rnw:166-167
###################################################
TR2 <- getMCMCtrees(2)


###################################################
### code chunk number 11: Running_coalescentMCMC.Rnw:172-173
###################################################
getMCMCstats()


###################################################
### code chunk number 12: Running_coalescentMCMC.Rnw:193-197
###################################################
tr <- TR[201:300]
THETA <- out[301:400, 2]
logLik0 <- mapply(dcoal, phy = tr, theta = THETA, log = TRUE)
summary(logLik0)


###################################################
### code chunk number 13: Running_coalescentMCMC.Rnw:201-206
###################################################
tr2 <- TR2[201:300]
THETA0 <- out2[301:400, 2]
RHO <- out2[301:400, 3]
logLik1 <- mapply(dcoal.time, phy = tr2, theta = THETA0, rho = RHO, log = TRUE)
summary(logLik1)


###################################################
### code chunk number 14: Running_coalescentMCMC.Rnw:210-212
###################################################
library(lattice)
print(histogram(~c(logLik0, logLik1) | gl(2, 100, labels = c("H0", "H1"))))


###################################################
### code chunk number 15: Running_coalescentMCMC.Rnw:218-220
###################################################
(MLE0 <- colMeans(out[301:400, 2, drop = FALSE]))
(MLE1 <- colMeans(out2[301:400, 2:3]))


###################################################
### code chunk number 16: Running_coalescentMCMC.Rnw:224-227
###################################################
D2.0 <- -2 * mean(mapply(dcoal, phy = tr, theta = MLE0, log = TRUE)) + 2
D2.1 <- -2 * mean(mapply(dcoal.time, phy = tr2, theta = MLE1[1],
                         rho = MLE1[2], log = TRUE)) + 4


###################################################
### code chunk number 17: Running_coalescentMCMC.Rnw:231-233
###################################################
mean(-2 * logLik0 + 2) - D2.0
mean(-2 * logLik1 + 4) - D2.1


###################################################
### code chunk number 18: Running_coalescentMCMC.Rnw:239-242
###################################################
x <- seq(0, 0.01, 0.0001)
y <- MLE1["theta0"] * exp(MLE1["rho"] * x)
plot(-x, y, "l", xlab = "Time", ylab = expression(Theta))


