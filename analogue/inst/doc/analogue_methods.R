### R code from vignette source 'analogue_methods.Rnw'

###################################################
### code chunk number 1: preliminary
###################################################
options("prompt" = "R> ", "continue" = "+ ")


###################################################
### code chunk number 2: loadPackage
###################################################
library("analogue")


###################################################
### code chunk number 3: loadData
###################################################
data(swapdiat, swappH, rlgh, package = "analogue")


###################################################
### code chunk number 4: joinData
###################################################
dat <- join(swapdiat, rlgh, verbose = TRUE)


###################################################
### code chunk number 5: convertToPercent
###################################################
swapdiat <- dat$swapdiat / 100
rlgh <- dat$rlgh / 100


###################################################
### code chunk number 6: applyMat
###################################################
swap.mat <- mat(swapdiat, swappH, method = "SQchord")


###################################################
### code chunk number 7: analogue_methods.Rnw:134-135
###################################################
swap.mat


###################################################
### code chunk number 8: summarySwapMat
###################################################
summary(swap.mat)


###################################################
### code chunk number 9: plot_mat
###################################################
opar <- par(mfrow = c(2,2))
plot(swap.mat)
par(opar)


###################################################
### code chunk number 10: printSwapMat
###################################################
getK(swap.mat)


###################################################
### code chunk number 11: analogue_methods.Rnw:168-169
###################################################
opar <- par(mfrow = c(2,2))
plot(swap.mat)
par(opar)


###################################################
### code chunk number 12: rlghPred
###################################################
rlgh.mat <- predict(swap.mat, rlgh, k = 10)
rlgh.mat


###################################################
### code chunk number 13: plot_recon
###################################################
reconPlot(rlgh.mat, use.labels = TRUE, ylab = "pH", xlab = "Depth (cm.)")


###################################################
### code chunk number 14: analogue_methods.Rnw:193-194
###################################################
reconPlot(rlgh.mat, use.labels = TRUE, ylab = "pH", xlab = "Depth (cm.)")


###################################################
### code chunk number 15: minDijRLGH
###################################################
rlgh.mdc <- minDC(rlgh.mat)


###################################################
### code chunk number 16: plot_minDC
###################################################
plot(rlgh.mdc, use.labels = TRUE, xlab = "Depth (cm.)")


###################################################
### code chunk number 17: analogue_methods.Rnw:215-216
###################################################
plot(rlgh.mdc, use.labels = TRUE, xlab = "Depth (cm.)")


###################################################
### code chunk number 18: swapAnalogue
###################################################
rlgh.ref <- rlgh[25:37, ]
swap.ana <- analog(swapdiat, rlgh.ref, method = "chord")
swap.ana


###################################################
### code chunk number 19: cma
###################################################
swap.cma <- cma(swap.ana)
swap.cma


###################################################
### code chunk number 20: analogue_methods.Rnw:245-246
###################################################
cma(swap.ana, cutoff = 0.5)


###################################################
### code chunk number 21: plot_cma
###################################################
plot(swap.cma)


###################################################
### code chunk number 22: analogue_methods.Rnw:260-261
###################################################
plot(swap.cma)


###################################################
### code chunk number 23: plot_dissim
###################################################
plot(dissim(swap.ana))


###################################################
### code chunk number 24: analogue_methods.Rnw:277-278
###################################################
plot(dissim(swap.ana))


###################################################
### code chunk number 25: mcarlo
###################################################
swap.mc <- mcarlo(swap.ana)
swap.mc


###################################################
### code chunk number 26: bootstrap
###################################################
set.seed(1234)
swap.boot <- bootstrap(swap.mat, n.boot = 100)
swap.boot


###################################################
### code chunk number 27: rmsep
###################################################
RMSEP(swap.boot, type = "standard")


###################################################
### code chunk number 28: analogue_methods.Rnw:329-331
###################################################
getK(swap.boot)
setK(swap.mat) <- getK(swap.boot)


###################################################
### code chunk number 29: clust_sites
###################################################
clust <- hclust(as.dist(swap.mat$Dij), method = "ward") #$
grps <- cutree(clust, k = 12)


###################################################
### code chunk number 30: roc
###################################################
swap.roc <- roc(swap.mat, groups = grps)
swap.roc


###################################################
### code chunk number 31: plot_roc
###################################################
opar <- par(mfrow = c(2,2))
plot(swap.roc)
par(opar)


###################################################
### code chunk number 32: analogue_methods.Rnw:370-371
###################################################
opar <- par(mfrow = c(2,2))
plot(swap.roc)
par(opar)


###################################################
### code chunk number 33: analogue_methods.Rnw:439-441
###################################################
dists1 <- distance(swapdiat, method = "bray")
dists2 <- distance(swapdiat, rlgh, method = "bray")


###################################################
### code chunk number 34: sample_specific
###################################################
set.seed(1234)
rlgh.boot <- predict(swap.mat, rlgh, bootstrap = TRUE, n.boot = 100)
reconPlot(rlgh.boot, use.labels = TRUE, ylab = "pH", xlab = "Depth (cm.)",
          display.error = "bars", predictions = "bootstrap")


###################################################
### code chunk number 35: analogue_methods.Rnw:468-470
###################################################
reconPlot(rlgh.boot, use.labels = TRUE, ylab = "pH", xlab = "Depth (cm.)",
          display.error = "bars", predictions = "bootstrap")


###################################################
### code chunk number 36: analogue_methods.Rnw:480-486
###################################################
set.seed(1234)
want <- sample(1:nrow(swapdiat), 67, replace = FALSE)
train <- swapdiat[-want, ]
train.env <- swappH[-want]
test <- swapdiat[want, ]
test.env <- swappH[want]


###################################################
### code chunk number 37: testset_boot
###################################################
train.mat <- mat(train, train.env, method = "SQchord")
test.boot <- bootstrap(train.mat, newdata = test,
                       newenv = test.env, n.boot = 100)
test.boot


###################################################
### code chunk number 38: analogue_methods.Rnw:504-510
###################################################
set.seed(9876)
want <- sample(nrow(test), 40)
opti <- test[-want, ]
opti.env <- test.env[-want]
test <- test[want, ]
test.env <- test.env[want]


###################################################
### code chunk number 39: opti
###################################################
opti.boot <- bootstrap(train.mat, newdata = opti, newenv = opti.env, n.boot = 100)
opti.boot


###################################################
### code chunk number 40: analogue_methods.Rnw:523-526
###################################################
use.k <- getK(opti.boot, prediction = TRUE, which = "model")
test.boot <- bootstrap(train.mat, newdata = test, newenv = test.env, k = use.k, n.boot = 100)
test.boot


###################################################
### code chunk number 41: subsetting
###################################################
dat <- join(swapdiat, rlgh, split = FALSE)
max.abb <- apply(dat, 2, max)
n.occ <- colSums(dat > 0)
spp.want <- which(max.abb >= 0.02 & n.occ >= 5)
swapdiat2 <- swapdiat[, spp.want]
rlgh2 <- rlgh[, spp.want]


###################################################
### code chunk number 42: mat_inferred
###################################################
plot(swap.mat, which = 1)


###################################################
### code chunk number 43: analogue_methods.Rnw:563-564
###################################################
plot(swap.mat, which = 1)


###################################################
### code chunk number 44: mat_resid
###################################################
plot(swap.mat, which = 2)


###################################################
### code chunk number 45: analogue_methods.Rnw:581-582
###################################################
plot(swap.mat, which = 2)


###################################################
### code chunk number 46: scree
###################################################
screeplot(swap.boot)


###################################################
### code chunk number 47: analogue_methods.Rnw:599-600
###################################################
screeplot(swap.boot)


