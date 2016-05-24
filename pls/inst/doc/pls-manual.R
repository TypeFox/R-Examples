### R code from vignette source 'pls-manual.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: pls-manual.Rnw:69-71
###################################################
pdf.options(pointsize=10)
options(digits = 4)


###################################################
### code chunk number 2: pls-manual.Rnw:306-307
###################################################
library(pls)


###################################################
### code chunk number 3: pls-manual.Rnw:334-337
###################################################
data(yarn)
data(oliveoil)
data(gasoline)


###################################################
### code chunk number 4: pls-manual.Rnw:345-351
###################################################
par(mar = c(2, 4, 0, 1) + 0.1)
matplot(t(gasoline$NIR), type = "l", lty = 1, ylab = "log(1/R)", xaxt = "n")
ind <- pretty(seq(from = 900, to = 1700, by = 2))
ind <- ind[ind >= 900 & ind <= 1700]
ind <- (ind - 898) / 2
axis(1, ind, colnames(gasoline$NIR)[ind])


###################################################
### code chunk number 5: pls-manual.Rnw:360-362
###################################################
gasTrain <- gasoline[1:50,]
gasTest <- gasoline[51:60,]


###################################################
### code chunk number 6: pls-manual.Rnw:365-366
###################################################
gas1 <- plsr(octane ~ NIR, ncomp = 10, data = gasTrain, validation = "LOO")


###################################################
### code chunk number 7: pls-manual.Rnw:371-372
###################################################
summary(gas1)


###################################################
### code chunk number 8: pls-manual.Rnw:381-382 (eval = FALSE)
###################################################
## plot(RMSEP(gas1), legendpos = "topright")


###################################################
### code chunk number 9: pls-manual.Rnw:387-389
###################################################
par(mar = c(4, 4, 2.5, 1) + 0.1)
plot(RMSEP(gas1), legendpos = "topright")


###################################################
### code chunk number 10: pls-manual.Rnw:408-409 (eval = FALSE)
###################################################
## plot(gas1, ncomp = 2, asp = 1, line = TRUE)


###################################################
### code chunk number 11: pls-manual.Rnw:415-417
###################################################
par(mar = c(4, 4, 2.5, 1) + 0.1)
plot(gas1, ncomp = 2, asp = 1, line = TRUE)


###################################################
### code chunk number 12: pls-manual.Rnw:431-432 (eval = FALSE)
###################################################
## plot(gas1, plottype = "scores", comps = 1:3)


###################################################
### code chunk number 13: pls-manual.Rnw:437-438
###################################################
plot(gas1, plottype = "scores", comps = 1:3)


###################################################
### code chunk number 14: pls-manual.Rnw:446-450
###################################################
par(mar = c(4, 4, 0.3, 1) + 0.1)
plot(gas1, "loadings", comps = 1:2, legendpos = "topleft",
     labels = "numbers", xlab = "nm")
abline(h = 0)


###################################################
### code chunk number 15: pls-manual.Rnw:464-465
###################################################
explvar(gas1)


###################################################
### code chunk number 16: pls-manual.Rnw:471-474 (eval = FALSE)
###################################################
## plot(gas1, "loadings", comps = 1:2, legendpos = "topleft",
##      labels = "numbers", xlab = "nm")
## abline(h = 0)


###################################################
### code chunk number 17: pls-manual.Rnw:483-484
###################################################
predict(gas1, ncomp = 2, newdata = gasTest)


###################################################
### code chunk number 18: pls-manual.Rnw:488-489
###################################################
RMSEP(gas1, newdata = gasTest)


###################################################
### code chunk number 19: pls-manual.Rnw:589-590
###################################################
dens1 <- plsr(density ~ NIR, ncomp = 5, data = yarn)


###################################################
### code chunk number 20: pls-manual.Rnw:594-596
###################################################
dim(oliveoil$sensory)
plsr(sensory ~ chemical, data = oliveoil)


###################################################
### code chunk number 21: pls-manual.Rnw:617-619
###################################################
trainind <- which(yarn$train == TRUE)
dens2 <- update(dens1, subset = trainind)


###################################################
### code chunk number 22: pls-manual.Rnw:623-624
###################################################
dens3 <- update(dens1, ncomp = 10)


###################################################
### code chunk number 23: pls-manual.Rnw:654-655
###################################################
olive1 <- plsr(sensory ~ chemical, scale = TRUE, data = oliveoil)


###################################################
### code chunk number 24: pls-manual.Rnw:660-661
###################################################
gas2 <- plsr(octane ~ msc(NIR), ncomp = 10, data = gasTrain)


###################################################
### code chunk number 25: pls-manual.Rnw:666-667 (eval = FALSE)
###################################################
## predict(gas2, ncomp = 3, newdata = gasTest)


###################################################
### code chunk number 26: pls-manual.Rnw:716-719
###################################################
gas2.cv <- crossval(gas2, segments = 10)
plot(MSEP(gas2.cv), legendpos="topright")
summary(gas2.cv, what = "validation")


###################################################
### code chunk number 27: pls-manual.Rnw:774-776 (eval = FALSE)
###################################################
## plot(gas1, plottype = "coef", ncomp=1:3, legendpos = "bottomleft",
##      labels = "numbers", xlab = "nm")


###################################################
### code chunk number 28: pls-manual.Rnw:780-783
###################################################
par(mar = c(4, 4, 2.5, 1) + 0.1)
plot(gas1, plottype = "coef", ncomp=1:3, legendpos = "bottomleft",
     labels = "numbers", xlab = "nm")


###################################################
### code chunk number 29: pls-manual.Rnw:815-817
###################################################
par(mar = c(4, 4, 0, 1) + 0.1)
plot(gas1, plottype = "correlation")


###################################################
### code chunk number 30: pls-manual.Rnw:886-887
###################################################
predict(gas1, ncomp = 2:3, newdata = gasTest[1:5,])


###################################################
### code chunk number 31: pls-manual.Rnw:899-900
###################################################
predict(gas1, comps = 2, newdata = gasTest[1:5,])


###################################################
### code chunk number 32: pls-manual.Rnw:917-918
###################################################
drop(predict(gas1, ncomp = 2:3, newdata = gasTest[1:5,]))


###################################################
### code chunk number 33: pls-manual.Rnw:956-957 (eval = FALSE)
###################################################
## predplot(gas1, ncomp = 2, newdata = gasTest, asp = 1, line = TRUE)


###################################################
### code chunk number 34: pls-manual.Rnw:963-965
###################################################
par(mar = c(4, 4, 2.5, 1))
predplot(gas1, ncomp = 2, newdata = gasTest, asp = 1, line = TRUE)


###################################################
### code chunk number 35: pls-manual.Rnw:1005-1006
###################################################
pls.options()


###################################################
### code chunk number 36: pls-manual.Rnw:1012-1013
###################################################
pls.options(plsralg = "oscorespls")


###################################################
### code chunk number 37: pls-manual.Rnw:1166-1175
###################################################
X <- gasTrain$NIR
Y <- gasTrain$octane
ncomp <- 5
cvPreds <- matrix(nrow = nrow(X), ncol = ncomp)
for (i in 1:nrow(X)) {
    fit <- simpls.fit(X[-i,], Y[-i], ncomp = ncomp, stripped = TRUE)
    cvPreds[i,] <- (X[i,] - fit$Xmeans) %*% drop(fit$coefficients) +
        fit$Ymeans
}


###################################################
### code chunk number 38: pls-manual.Rnw:1178-1179
###################################################
sqrt(colMeans((cvPreds - Y)^2))


