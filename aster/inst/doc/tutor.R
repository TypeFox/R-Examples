### R code from vignette source 'tutor.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: setup
###################################################
library(aster)
packageDescription("aster")$Version
data(echinacea)


###################################################
### code chunk number 2: data-names
###################################################
names(echinacea)


###################################################
### code chunk number 3: data-classes
###################################################
sapply(echinacea, class)


###################################################
### code chunk number 4: reshape
###################################################
vars <- c("ld02", "ld03", "ld04", "fl02", "fl03", "fl04",
    "hdct02", "hdct03", "hdct04")
redata <- reshape(echinacea, varying = list(vars),
     direction = "long", timevar = "varb", times = as.factor(vars),
     v.names = "resp")
names(redata)


###################################################
### code chunk number 5: show-redata-varb
###################################################
class(redata$varb)
levels(redata$varb)


###################################################
### code chunk number 6: reshape-explain
###################################################
redata[42, ]


###################################################
### code chunk number 7: reshape-too
###################################################
redata <- data.frame(redata, root = 1)
names(redata)


###################################################
### code chunk number 8: graph
###################################################
pred <- c(0, 1, 2, 1, 2, 3, 4, 5, 6)
fam <- c(1, 1, 1, 1, 1, 1, 3, 3, 3)
sapply(fam.default(), as.character)[fam]


###################################################
### code chunk number 9: fit-1
###################################################
aout1 <- aster(resp ~ varb + nsloc + ewloc + pop,
    pred, fam, varb, id, root, data = redata)
summary(aout1, show.graph = TRUE)


###################################################
### code chunk number 10: signif.stars
###################################################
options(show.signif.stars = FALSE)


###################################################
### code chunk number 11: fit-2
###################################################
aout2 <- aster(resp ~ varb + nsloc + ewloc,
    pred, fam, varb, id, root, data = redata)


###################################################
### code chunk number 12: fit-2-anova
###################################################
anova(aout2, aout1)


###################################################
### code chunk number 13: make-hdct
###################################################
hdct <- grep("hdct", as.character(redata$varb))
hdct <- is.element(seq(along = redata$varb), hdct)
redata <- data.frame(redata, hdct = as.integer(hdct))
names(redata)


###################################################
### code chunk number 14: fit-3
###################################################
aout3 <- aster(resp ~ varb + nsloc + ewloc + pop * hdct,
    pred, fam, varb, id, root, data = redata)
summary(aout3)


###################################################
### code chunk number 15: fit-4
###################################################
aout4 <- aster(resp ~ varb + nsloc + ewloc + pop * hdct - pop,
    pred, fam, varb, id, root, data = redata)
summary(aout4)


###################################################
### code chunk number 16: fit-4-anova
###################################################
anova(aout2, aout4, aout3)


###################################################
### code chunk number 17: conf-level
###################################################
conf.level <- 0.95
crit <- qnorm((1 + conf.level) / 2)


###################################################
### code chunk number 18: beta-4-one
###################################################
fred <- summary(aout4)
dimnames(fred$coef)
fred$coef["popEriley:hdct", "Estimate"] +
c(-1, 1) * crit * fred$coef["popEriley:hdct", "Std. Error"]


###################################################
### code chunk number 19: beta-4-one-fish
###################################################
aout4$coef[12] + c(-1, 1) * crit * sqrt(solve(aout4$fish)[12, 12])


###################################################
### code chunk number 20: beta-4-two-fish
###################################################
inv.fish.info <- solve(aout4$fish)
(aout4$coef[12] - aout4$coef[13]) + c(-1, 1) * crit *
sqrt(inv.fish.info[12, 12] + inv.fish.info[13, 13] - 2 * inv.fish.info[12, 13])


###################################################
### code chunk number 21: predict-newdata
###################################################
newdata <- data.frame(pop = levels(echinacea$pop))
for (v in vars)
    newdata[[v]] <- 1
newdata$root <- 1
newdata$ewloc <- 0
newdata$nsloc <- 0


###################################################
### code chunk number 22: predict-newdata-reshape
###################################################
renewdata <- reshape(newdata, varying = list(vars),
     direction = "long", timevar = "varb", times = as.factor(vars),
     v.names = "resp")
hdct <- grep("hdct", as.character(renewdata$varb))
hdct <- is.element(seq(along = renewdata$varb), hdct)
renewdata <- data.frame(renewdata, hdct = as.integer(hdct))
names(redata)
names(renewdata)


###################################################
### code chunk number 23: make-amat
###################################################
nind <- nrow(newdata)
nnode <- length(vars)
amat <- array(0, c(nind, nnode, nind))
for (i in 1:nind)
    amat[i , grep("hdct", vars), i] <- 1


###################################################
### code chunk number 24: varphi-4-red-fish
###################################################
foo <- predict(aout4, varvar = varb, idvar = id, root = root,
    newdata = renewdata, se.fit = TRUE, amat = amat,
    parm.type = "canon")
bar <- cbind(foo$fit, foo$se.fit)
dimnames(bar) <- list(as.character(newdata$pop), c("Estimate", "Std. Error"))
print(bar)


###################################################
### code chunk number 25: theta-4-blue-fish
###################################################
foo <- predict(aout4, varvar = varb, idvar = id, root = root,
    newdata = renewdata, se.fit = TRUE, amat = amat,
    parm.type = "canon", model.type = "cond")
bar <- cbind(foo$fit, foo$se.fit)
dimnames(bar) <- list(as.character(newdata$pop), c("Estimate", "Std. Error"))
print(bar)


###################################################
### code chunk number 26: make-bmat
###################################################
bmat <- array(0, c(nind, nnode, nind))
for (i in 1:nind)
    bmat[i , grep("ld", vars), i] <- 1


###################################################
### code chunk number 27: bmat-4-fish-knife
###################################################
foo <- predict(aout4, varvar = varb, idvar = id, root = root,
    newdata = renewdata, se.fit = TRUE, amat = bmat,
    parm.type = "canon")
bar <- cbind(foo$fit, foo$se.fit)
dimnames(bar) <- list(as.character(newdata$pop), c("Estimate", "Std. Error"))
print(bar)
foo <- predict(aout4, varvar = varb, idvar = id, root = root,
    newdata = renewdata, se.fit = TRUE, amat = bmat,
    parm.type = "canon", model.type = "cond")
bar <- cbind(foo$fit, foo$se.fit)
dimnames(bar) <- list(as.character(newdata$pop), c("Estimate", "Std. Error"))
print(bar)


###################################################
### code chunk number 28: bmat-4-fish-fash
###################################################
foo <- predict(aout3, varvar = varb, idvar = id, root = root,
    newdata = renewdata, se.fit = TRUE, amat = bmat,
    parm.type = "canon")
bar <- cbind(foo$fit, foo$se.fit)
dimnames(bar) <- list(as.character(newdata$pop), c("Estimate", "Std. Error"))
print(bar)


###################################################
### code chunk number 29: tau-4-amat
###################################################
pout3 <- predict(aout3, varvar = varb, idvar = id, root = root,
    newdata = renewdata, se.fit = TRUE, amat = amat)
pout4 <- predict(aout4, varvar = varb, idvar = id, root = root,
    newdata = renewdata, se.fit = TRUE, amat = amat)


###################################################
### code chunk number 30: fig1plot
###################################################
popnames <- as.character(newdata$pop)
fit3 <- pout3$fit
fit4 <- pout4$fit
i <- seq(along = popnames)
foo <- 0.1
y4top <- fit4 + crit * pout4$se.fit
y4bot <- fit4 - crit * pout4$se.fit
y3top <- fit3 + crit * pout3$se.fit
y3bot <- fit3 - crit * pout3$se.fit
plot(c(i - 1.5 * foo, i - 1.5 * foo, i + 1.5 * foo, i + 1.5 * foo),
    c(y4top, y4bot, y3top, y3bot), type = "n", axes = FALSE,
    xlab = "", ylab = "")
segments(i - 1.5 * foo, y4bot, i - 1.5 * foo, y4top)
segments(i - 2.5 * foo, y4bot, i - 0.5 * foo, y4bot)
segments(i - 2.5 * foo, y4top, i - 0.5 * foo, y4top)
segments(i - 2.5 * foo, fit4, i - 0.5 * foo, fit4)
segments(i + 1.5 * foo, y3bot, i + 1.5 * foo, y3top, lty = 2)
segments(i + 2.5 * foo, y3bot, i + 0.5 * foo, y3bot)
segments(i + 2.5 * foo, y3top, i + 0.5 * foo, y3top)
segments(i + 2.5 * foo, fit3, i + 0.5 * foo, fit3)
axis(side = 2)
title(ylab = "unconditional mean value parameter")
axis(side = 1, at = i, labels = popnames)
title(xlab = "population")


###################################################
### code chunk number 31: fig1
###################################################
popnames <- as.character(newdata$pop)
fit3 <- pout3$fit
fit4 <- pout4$fit
i <- seq(along = popnames)
foo <- 0.1
y4top <- fit4 + crit * pout4$se.fit
y4bot <- fit4 - crit * pout4$se.fit
y3top <- fit3 + crit * pout3$se.fit
y3bot <- fit3 - crit * pout3$se.fit
plot(c(i - 1.5 * foo, i - 1.5 * foo, i + 1.5 * foo, i + 1.5 * foo),
    c(y4top, y4bot, y3top, y3bot), type = "n", axes = FALSE,
    xlab = "", ylab = "")
segments(i - 1.5 * foo, y4bot, i - 1.5 * foo, y4top)
segments(i - 2.5 * foo, y4bot, i - 0.5 * foo, y4bot)
segments(i - 2.5 * foo, y4top, i - 0.5 * foo, y4top)
segments(i - 2.5 * foo, fit4, i - 0.5 * foo, fit4)
segments(i + 1.5 * foo, y3bot, i + 1.5 * foo, y3top, lty = 2)
segments(i + 2.5 * foo, y3bot, i + 0.5 * foo, y3bot)
segments(i + 2.5 * foo, y3top, i + 0.5 * foo, y3top)
segments(i + 2.5 * foo, fit3, i + 0.5 * foo, fit3)
axis(side = 2)
title(ylab = "unconditional mean value parameter")
axis(side = 1, at = i, labels = popnames)
title(xlab = "population")


###################################################
### code chunk number 32: make-cout
###################################################
cout4 <- aster(resp ~ varb + nsloc + ewloc + pop * hdct - pop,
    pred, fam, varb, id, root, data = redata, type = "cond")
pcout4 <- predict(cout4, varvar = varb, idvar = id, root = root,
    newdata = renewdata, se.fit = TRUE, amat = amat)


###################################################
### code chunk number 33: fig2plot
###################################################
popnames <- as.character(newdata$pop)
fit3 <- pcout4$fit
fit4 <- pout4$fit
i <- seq(along = popnames)
foo <- 0.1
y4top <- fit4 + crit * pout4$se.fit
y4bot <- fit4 - crit * pout4$se.fit
y3top <- fit3 + crit * pcout4$se.fit
y3bot <- fit3 - crit * pcout4$se.fit
plot(c(i - 1.5 * foo, i - 1.5 * foo, i + 1.5 * foo, i + 1.5 * foo),
    c(y4top, y4bot, y3top, y3bot), type = "n", axes = FALSE,
    xlab = "", ylab = "")
segments(i - 1.5 * foo, y4bot, i - 1.5 * foo, y4top)
segments(i - 2.5 * foo, y4bot, i - 0.5 * foo, y4bot)
segments(i - 2.5 * foo, y4top, i - 0.5 * foo, y4top)
segments(i - 2.5 * foo, fit4, i - 0.5 * foo, fit4)
segments(i + 1.5 * foo, y3bot, i + 1.5 * foo, y3top, lty = 2)
segments(i + 2.5 * foo, y3bot, i + 0.5 * foo, y3bot)
segments(i + 2.5 * foo, y3top, i + 0.5 * foo, y3top)
segments(i + 2.5 * foo, fit3, i + 0.5 * foo, fit3)
axis(side = 2)
title(ylab = "unconditional mean value parameter")
axis(side = 1, at = i, labels = popnames)
title(xlab = "population")


###################################################
### code chunk number 34: fig2
###################################################
popnames <- as.character(newdata$pop)
fit3 <- pcout4$fit
fit4 <- pout4$fit
i <- seq(along = popnames)
foo <- 0.1
y4top <- fit4 + crit * pout4$se.fit
y4bot <- fit4 - crit * pout4$se.fit
y3top <- fit3 + crit * pcout4$se.fit
y3bot <- fit3 - crit * pcout4$se.fit
plot(c(i - 1.5 * foo, i - 1.5 * foo, i + 1.5 * foo, i + 1.5 * foo),
    c(y4top, y4bot, y3top, y3bot), type = "n", axes = FALSE,
    xlab = "", ylab = "")
segments(i - 1.5 * foo, y4bot, i - 1.5 * foo, y4top)
segments(i - 2.5 * foo, y4bot, i - 0.5 * foo, y4bot)
segments(i - 2.5 * foo, y4top, i - 0.5 * foo, y4top)
segments(i - 2.5 * foo, fit4, i - 0.5 * foo, fit4)
segments(i + 1.5 * foo, y3bot, i + 1.5 * foo, y3top, lty = 2)
segments(i + 2.5 * foo, y3bot, i + 0.5 * foo, y3bot)
segments(i + 2.5 * foo, y3top, i + 0.5 * foo, y3top)
segments(i + 2.5 * foo, fit3, i + 0.5 * foo, fit3)
axis(side = 2)
title(ylab = "unconditional mean value parameter")
axis(side = 1, at = i, labels = popnames)
title(xlab = "population")


###################################################
### code chunk number 35: make-theta
###################################################
theta.hat <- predict(aout4, model.type = "cond", parm.type = "canon")
theta.hat <- matrix(theta.hat, nrow = nrow(aout4$x), ncol = ncol(aout4$x))
fit.hat <- pout4$fit
beta.hat <- aout4$coefficients


###################################################
### code chunk number 36: make-root
###################################################
root <- aout4$root
modmat <- aout4$modmat
modmat.pred <- pout4$modmat
x.pred <- matrix(1, nrow = dim(modmat.pred)[1], ncol = dim(modmat.pred)[2])
root.pred <- x.pred


###################################################
### code chunk number 37: doit
###################################################
set.seed(42)
nboot <- 100
cover <- matrix(0, nboot, length(fit.hat))
for (iboot in 1:nboot) {
    xstar <- raster(theta.hat, pred, fam, root)
    aout4star <- aster(xstar, root, pred, fam, modmat, beta.hat)
    pout4star <- predict(aout4star, x.pred, root.pred, modmat.pred,
        amat, se.fit = TRUE)
    upper <- pout4star$fit + crit * pout4star$se.fit
    lower <- pout4star$fit - crit * pout4star$se.fit
    cover[iboot, ] <- as.numeric(lower <= fit.hat & fit.hat <= upper)
}
pboot <- apply(cover, 2, mean)
pboot.se <- sqrt(pboot * (1 - pboot) / nboot)
cbind(pboot, pboot.se)


