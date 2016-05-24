### R code from vignette source 'twoACexamples.Rnw'

###################################################
### code chunk number 1: Initialize
###################################################

RUN <- FALSE    #redo computations and write .RData files
## Change options:
op <- options() ## To be able to reset settings
options("digits" = 7)
options("htmlhelp" = TRUE)
## options("width" = 75)
options("SweaveHooks" = list(fig=function()
        par(mar=c(4,4,1.5,1)+.5)))
options(continue=" ")



###################################################
### code chunk number 2: example1-1
###################################################
library(sensR)
fit <- twoAC(c(2,2,6))
fit


###################################################
### code chunk number 3: example1-2
###################################################
n <- c(2, 2, 6)
gamma <- cumsum(n / sum(n))
z <- qnorm(gamma)[-3]
z <- z * sqrt(2)
(tau <- (z[2] - z[1])/2)
(d <- -z[1] - tau)


###################################################
### code chunk number 4: example2
###################################################
twoAC(c(2, 2, 6), d.prime0 = 0, conf.level = 0.95,
      statistic = "likelihood", alternative = "greater")


###################################################
### code chunk number 5: example2-fig
###################################################
getOption("SweaveHooks")[["fig"]]()
pr <- profile(fit)
plot(pr)
z <- pr$d.prime$d.prime
w <- (coef(fit)[2,1] - z) / coef(fit)[2,2]
lines(z, exp(-w^2/2), lty = 2)
## pl <- plot(pr, fig = FALSE)
## abline(h = attr(pr, "limits"))
## text(2.7, .17, "95% limit")
## text(2.7, .06, "99% limit")



###################################################
### code chunk number 6: power-example
###################################################
twoACpwr(tau = 0.5, d.prime = 1, size = 20, d.prime0 = 0, alpha = 0.05,
         alternative = "two.sided", tol = 1e-5)


###################################################
### code chunk number 7: power-example2
###################################################
twoACpwr(tau = 0.5, d.prime = 1, size = 20, d.prime0 = 0, alpha = 0.05,
         alternative = "two.sided", tol = 0)


###################################################
### code chunk number 8: example2-1
###################################################
library(ordinal)
response <- gl(3,1)
fit.clm <- clm(response ~ 1, weights = c(2, 2, 6), link = "probit")
(tab <- coef(summary(fit.clm)))


###################################################
### code chunk number 9: example2-2
###################################################
theta <- tab[,1]
(tau <- (theta[2] - theta[1])/sqrt(2))
(d.prime <- (-theta[2] - theta[1])/sqrt(2))


###################################################
### code chunk number 10: example2-3
###################################################
## SE <- sqrt(diag(vcov(fit.clm)))
VCOV <- vcov(fit.clm)
(se.tau <- sqrt((VCOV[1,1] + VCOV[2,2] - 2*VCOV[2,1])/2))
(se.d.prime <- sqrt((VCOV[1,1] + VCOV[2,2] + 2*VCOV[2,1])/2))
##
## (SE[2] + SE[1]) * sqrt(2) / 2
## theta <- fit.clm$Theta
## ## Delta <- theta[2] - theta[1]
## ## theta0 <- sum(theta)/2
## ## (d <- -sum(fit.clm$Theta * sqrt(2))/2)
## ## theta[2] - theta0
## ## diff(fit.clm$Theta * sqrt(2))/2
## ###
## (tau <- (theta[2] - theta[1])/sqrt(2))
## (d <- (-theta[2] - theta[1])/sqrt(2))


###################################################
### code chunk number 11: example2-4
###################################################
clm2twoAC(fit.clm)


###################################################
### code chunk number 12: example3-1
###################################################
n.women <- c(2, 2, 6)*10
## n.women <- c(2, 2, 6)*10
n.men <- c(1, 2, 7)*10
## n.men <- c(2, 2, 6)*10
wt <- c(n.women, n.men)
response <- gl(3,1, length = 6)
gender <- gl(2, 3, labels = c("women", "men"))
fm2 <- clm(response ~ gender, weights = wt, link = "probit")
(tab2 <- coef(summary(fm2)))


###################################################
### code chunk number 13: example3-2
###################################################
theta <- fm2$alpha
(tau <- (theta[2] - theta[1])/sqrt(2))
(d.women <- (-theta[2] - theta[1])/sqrt(2))
(d.men <- d.women + fm2$beta * sqrt(2))


###################################################
### code chunk number 14: example3-2_2
###################################################
clm2twoAC(fm2)


###################################################
### code chunk number 15: example3-3
###################################################
fm3 <- update(fm2,~.-gender)
anova(fm2, fm3)


###################################################
### code chunk number 16: example3-4
###################################################
confint(fm2) *  sqrt(2)
## pr <- profile(fm2, alpha = 1e-5)
## plot(pr, n=1e3)


###################################################
### code chunk number 17: example3-5
###################################################
logLik(fm2)
tw <- twoAC(n.women)
tm <- twoAC(n.men)
(LR <- 2*(tw$logLik + tm$logLik - fm2$logLik))
pchisq(LR, 1, lower.tail=FALSE)


###################################################
### code chunk number 18: example3-7
###################################################
freq <- matrix(fitted(fm2), nrow=2, byrow=TRUE) * 100
Obs <- matrix(wt, nrow=2, byrow=TRUE)
(X2 <- sum((Obs - freq)^2 / freq))
pchisq(X2, df=1, lower.tail=FALSE)


###################################################
### code chunk number 19: readData
###################################################
## load("C:/Users/rhbc/Documents/Rpackages/sensR/pkg/inst/doc/HSdata.RData")
load("HSdata.RData")
repData <- hs.long


###################################################
### code chunk number 20: Example6-1
###################################################
fm3.agq <- clmm2(preference ~ reference, random = cons, nAGQ = 10,
                data = repData, link = "probit", Hess = TRUE)
summary(fm3.agq)


###################################################
### code chunk number 21: Example6-2
###################################################
clm2twoAC(fm3.agq)


###################################################
### code chunk number 22: example6-3
###################################################
fm3.agq$stDev * sqrt(2)


###################################################
### code chunk number 23: example6-4 (eval = FALSE)
###################################################
## pr <- profile(fm3.agq, range = c(.7, 1.8))
## ## save(pr, file = "profileSigma.RData")


###################################################
### code chunk number 24: example6-4b
###################################################
## load("C:/Users/rhbc/Documents/Rpackages/sensR/pkg/inst/doc/profileSigma.RData")
load("profileSigma.RData")


###################################################
### code chunk number 25: example6-4-fig
###################################################
getOption("SweaveHooks")[["fig"]]()
plpr <- plot(pr, fig = FALSE)
plot(sqrt(2) * plpr$stDev$x, plpr$stDev$y, type = "l",
     xlab = expression(sigma[delta]),
     ylab = "Relative profile likelihood", xlim = c(1, 2.5),
     axes = FALSE)
axis(1); axis(2, las = 1)
abline(h = attr(plpr, "limits"))
text(2.4, .17, "95% limit")
text(2.4, .06, "99% limit")


###################################################
### code chunk number 26: example6-5
###################################################
confint(pr, level = 0.95) * sqrt(2)


###################################################
### code chunk number 27: example6-6
###################################################
newdat <- expand.grid(preference = factor(c("A", "N", "B"),
                        levels = c("A", "N", "B"), ordered = TRUE),
                      reference = factor(c("A", "B")))
pred <- predict(fm3.agq, newdata = newdat)


###################################################
### code chunk number 28: example6-7
###################################################
q95.refA <- diff(c(0, pnorm(fm3.agq$Theta - qnorm(1-.05) *
                            fm3.agq$stDev), 1))
q05.refA <- diff(c(0, pnorm(fm3.agq$Theta - qnorm(.05) *
                            fm3.agq$stDev), 1))
q95.refB <- diff(c(0, pnorm(fm3.agq$Theta - fm3.agq$beta -
                            qnorm(1-.05) * fm3.agq$stDev), 1))
q05.refB <- diff(c(0, pnorm(fm3.agq$Theta - fm3.agq$beta -
                            qnorm(.05) * fm3.agq$stDev), 1))


###################################################
### code chunk number 29: example6-8
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar = c(0, 2, 0, .5) + 0.5)
plot(1:3, pred[1:3], ylim = c(0, 1), axes = FALSE,
     xlab = "", ylab = "", pch = 19) ## ref A
axis(1, lwd.ticks = 0, at = c(1, 3), labels = c("", ""))
axis(2, las = 1)
points(1:3, pred[4:6], pch = 1) ## ref B
lines(1:3, pred[1:3]) ## ref A
lines(1:3, pred[4:6], lty = 2) ## ref B
text(2, 0.6, "Average consumer")
legend("topright", c("Reference A", "Reference B"),
       lty = 1:2, pch = c(19, 1), bty = "n")


###################################################
### code chunk number 30: example6-9
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar = c(0, 2, 0, .5) + 0.5)
plot(1:3, q05.refA, ylim = c(0, 1), axes = FALSE,
     xlab = "", ylab = "", pch = 19) ## ref A
## axis(1, at = 1:3, labels = c("prefer A", "no preference", "prefer B"))
axis(1, lwd.ticks = 0, at = c(1, 3), labels = c("", ""))
axis(2, las = 1)
points(1:3, q05.refB, pch = 1) ## ref B
lines(1:3, q05.refA) ## ref A
lines(1:3, q05.refB, lty = 2) ## ref B
text(2, 0.6, "5th percentile consumer")


###################################################
### code chunk number 31: example6-10
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar = c(2, 2, 0, .5) + 0.5)
plot(1:3, q95.refA, ylim = c(0, 1), axes = FALSE,
     xlab = "", ylab = "", pch = 19) ## ref A
axis(1, at = 1:3, labels = c("prefer A", "no preference", "prefer B"))
axis(2, las = 1)
points(1:3, q95.refB, pch = 1) ## ref B
lines(1:3, q95.refA) ## ref A
lines(1:3, q95.refB, lty = 2) ## ref B
text(2, 0.6, "95th percentile consumer")


###################################################
### code chunk number 32: sesstionInfo
###################################################
sessionInfo()


