### R code from vignette source 'methodology.Rnw'

###################################################
### code chunk number 1: Initialize
###################################################

## Load common packages, functions and set settings:
library(sensR)
##
RUN <- FALSE    #redo computations and write .RData files
## Change options:
op <- options() ## To be able to reset settings
options("digits" = 7)
options(help_type = "html")
options("width" = 75)
options("SweaveHooks" = list(fig=function()
        par(mar=c(4,4,.5,0)+.5)))
options(continue=" ")

G.test <- function(x, ...) {
    ct <- eval.parent(chisq.test(x, correct = FALSE))
    ct$data.name <- deparse(substitute(x))
    o <- ct$observed
    e <- ct$expected
    ct$statistic <- 2 * sum(ifelse(o == 0, 0, o * log(o/e)))
    names(ct$statistic) <- "G-squared"
    ct$p.value <- with(ct, pchisq(statistic, parameter, lower = FALSE))
    ct$method <- "Likelihood ratio Chi-squared test"
    return(ct)
}



###################################################
### code chunk number 2: scaleRelations
###################################################

pd <- seq(0, 1-1e-5, length = 1e2)
coef.twoAFC <- coef(rescale(pd = pd, method = "twoAFC"))
coef.threeAFC <- coef(rescale(pd = pd, method = "threeAFC"))
pd <- seq(0, 1-1e-2, length = 1e2)
coef.triangle <- coef(rescale(pd = pd, method = "triangle"))
coef.duotrio <- coef(rescale(pd = pd, method = "duotrio"))

## OK:
tail(coef.duotrio)
tail(coef.twoAFC)
tail(coef.threeAFC)
tail(coef.triangle)

dPrime <- seq(1e-5, 6, len = 1e2)
gd1 <- psyderiv(dPrime, method = "twoAFC")
gd2 <- psyderiv(dPrime, method = "threeAFC")
gd3 <- psyderiv(dPrime, method = "duotrio")
gd4 <- psyderiv(dPrime, method = "triangle")



###################################################
### code chunk number 3: psyFunsFig
###################################################
getOption("SweaveHooks")[["fig"]]()

plot(c(0, 6), c(.3, 1), type = "n", xlab = "d-prime",
     ylab = "P(correct answer)", axes = FALSE)
axis(1)
axis(2, las = 1)
with(coef.twoAFC, lines(d.prime, pc, lty = 1))
with(coef.threeAFC, lines(d.prime, pc, lty = 2))
with(coef.duotrio, lines(d.prime, pc, lty = 3))
with(coef.triangle, lines(d.prime, pc, lty = 4))
legend("topleft", legend = c("2-AFC", "3-AFC", "duo-trio",
                        "triangle"), lty = 1:4, bty = "n")



###################################################
### code chunk number 4: psyFunsFig_dprimePd
###################################################
getOption("SweaveHooks")[["fig"]]()

plot(c(0, 6), c(0, 1), type = "n", xlab = "d-prime",
     ylab = "P(discrimination)", axes = FALSE)
axis(1)
axis(2, las = 1)
with(coef.twoAFC, lines(d.prime, pd, lty = 1))
with(coef.threeAFC, lines(d.prime, pd, lty = 2))
with(coef.duotrio, lines(d.prime, pd, lty = 3))
with(coef.triangle, lines(d.prime, pd, lty = 4))
legend("topleft", legend = c("2-AFC", "3-AFC", "duo-trio",
                        "triangle"), lty = 1:4, bty = "n")



###################################################
### code chunk number 5: psyFunsFig_PcPd
###################################################
getOption("SweaveHooks")[["fig"]]()

plot(c(0.3, 1), c(0, 1), type = "n",
     ylab = "P(correct answer)",
     xlab = "P(discrimination)", axes = FALSE)
axis(1); axis(2, las = 1)
segments(c(1/3, .5), 0, 1, 1, lty = 1:2)
legend("topleft", legend = c("3-AFC and triangle", "2-AFC and duo-trio"),
       lty = 1:2, bty = "n")



###################################################
### code chunk number 6: psyFunFig_psyDeriv
###################################################
getOption("SweaveHooks")[["fig"]]()

plot(c(0, 6), c(0, .31), type = "n", axes = FALSE,
     xlab = "d-prime",
     ylab = "Grad(psychometic function)")
axis(1); axis(2, las = 1)
lines(dPrime, gd1, lty = 1)
lines(dPrime, gd2, lty = 2)
lines(dPrime, gd3, lty = 3)
lines(dPrime, gd4, lty = 4)
legend("topright", legend = c("2-AFC", "3-AFC", "duo-trio",
                     "triangle"), lty = 1:4, bty = "n")



###################################################
### code chunk number 7: rescaleExample
###################################################
rescale(pc = 0.25, method = "triangle")


###################################################
### code chunk number 8: rescaleExample2
###################################################
rescale(pd = 0.2, std.err = 0.12, method = "triangle")


###################################################
### code chunk number 9: discrimExample1
###################################################
discrim(10, 15, method = "threeAFC", statistic = "likelihood")


###################################################
### code chunk number 10: discrimExample2
###################################################
discrim(4, 15, method = "threeAFC", test = "similarity", pd0 = 0.2,
        statistic = "exact")


###################################################
### code chunk number 11: discrimExample3
###################################################
getOption("SweaveHooks")[["fig"]]()
fm1 <- discrim(10, 15, method = "threeAFC", statistic = "exact")
confint(fm1)
plot(profile(fm1))


###################################################
### code chunk number 12: discrimExample3Fig
###################################################
getOption("SweaveHooks")[["fig"]]()
fm1 <- discrim(10, 15, method = "threeAFC", statistic = "exact")
confint(fm1)
plot(profile(fm1))


###################################################
### code chunk number 13: criticalValueExample1
###################################################
1 - pbinom(q = 15 - 1, size = 20, prob = 0.5)
1 - pbinom(q = 14 - 1, size = 20, prob = 0.5)


###################################################
### code chunk number 14: criticalValueExample2
###################################################
i <- 0
while (1 - pbinom(q = i, size = 20, prob = 0.5) > 0.05)
  {
    i <- i + 1
  }
i + 1


###################################################
### code chunk number 15: criticalValueExample3
###################################################
findcr(sample.size = 20, alpha = 0.05, p0 = 0.5)


###################################################
### code chunk number 16: differencePowerExample1
###################################################
1 - pbinom(q = 15 - 1, size = 20, prob = 3/4)


###################################################
### code chunk number 17: differencePowerExample2
###################################################
discrimPwr(pdA = 0.5, sample.size = 20, alpha = 0.05,
           pGuess = 1/2)


###################################################
### code chunk number 18: differencePowerExample3
###################################################
discrimPwr(pdA = 0.5, pd0 = 0.1, sample.size = 20, alpha = 0.05,
           pGuess = 1/2)


###################################################
### code chunk number 19: similarityPowerExample1
###################################################
discrimPwr(pdA = 0, pd0 = 1/3, sample.size = 100, alpha = 0.05,
           pGuess = 1/2, test = "similarity")


###################################################
### code chunk number 20: similarityPowerExample1
###################################################
discrimPwr(pdA = 1/5, pd0 = 1/3, sample.size = 100, alpha = 0.05,
           pGuess = 1/2, test = "similarity")


###################################################
### code chunk number 21: powerSimul
###################################################
## Simulating power of a duo-trio no-difference test:
set.seed(12345)
xx <- rbinom(10000, 20, 3/4)
sum(xx == 0) ## 0
pp <- 1 - pbinom(xx-1, 20, 1/2)
sum(pp <= 0.05)
(power <- mean(pp <= 0.05)) # 0.6184

## Uncertainty (standard error of power):
(se.power <- sqrt(power * (1 - power) / 1e4))
## confidence interval:
power + c(-1,1) * qnorm(.975) * se.power


###################################################
### code chunk number 22: sampleSizeFigure
###################################################
ss <- 275:325
(pd <- coef(rescale(d.prime = .9, method = "triangle"))$pd)
pwr <- sapply(ss, function(x)
              discrimPwr(pdA = pd, sample.size = x, pGuess = 1/3))
pwrN <- sapply(ss, function(x)
               discrimPwr(pdA = pd, sample.size = x, pGuess = 1/3,
                          statistic = "normal"))



###################################################
### code chunk number 23: sampleSizeFigurePlot
###################################################
getOption("SweaveHooks")[["fig"]]()

plot(range(ss), c(0.74, max(pwrN)), type = "n", xlab = "sample size",
     ylab = "power", axes = FALSE)
lines(ss, pwr, type = "l")
points(ss, pwr, pch = 20, cex = .8)
abline(h = 0.8, lty = 2)
lines(ss, pwrN, col = "blue")
legend("topleft", legend = c("exact binomial", "normal approx.",
                    "target power"), lty = c(1, 1, 2),
       pch = c(20, NA, NA), col = c(1, "blue", 1), pt.cex = .8)
## axis(2, at = c(0.72, 0.75, 0.80, 0.85))
axis(2, las = 1)
axis(1, at = c(275, 297, 297+21, 325))
segments(297, 0.6, 297,
         discrimPwr(pd, sample.size = 297, pGuess = 1/3), lty = 3)
segments(297+21, 0.6, 297+21, discrimPwr(pd, sample.size = 297+21,
                                         pGuess = 1/3), lty = 3)



###################################################
### code chunk number 24: xcrFig
###################################################
getOption("SweaveHooks")[["fig"]]()
xcr <- sapply(ss, function(ss) findcr(ss, p0 = 1/3))
plot(ss, xcr, type = "l", axes = FALSE,
     ylab = "Critical value", xlab = "Sample size")
points(ss, xcr, pch = 20)
axis(1); axis(2)



###################################################
### code chunk number 25: alphaFig
###################################################
getOption("SweaveHooks")[["fig"]]()
aa <- sapply(ss, function(n)
             1 - pbinom(findcr(n, p0 = 1/3) - 1, n, 1/3))
plot(ss, aa, pch = 20, ylim = c(0.03, .05), axes = FALSE,
     ylab = "Actual alpha-level", xlab = "Sample size")
lines(ss, aa, lty = 1)
axis(1); axis(2)



###################################################
### code chunk number 26: exactSSexample
###################################################
(pd <- coef(rescale(d.prime = .9, method = "triangle"))$pd)
discrimSS(pdA = pd, pd0 = 0, target.power = 0.8, alpha = 0.05, pGuess
          = 1/3, test = "difference", statistic = "exact")


###################################################
### code chunk number 27: normalSSexample
###################################################
discrimSS(pdA = pd, pd0 = 0, target.power = 0.8, alpha = 0.05, pGuess
          = 1/3, test = "difference", statistic = "normal")


###################################################
### code chunk number 28: methodology.Rnw:1828-1831 (eval = FALSE)
###################################################
## 
## 
## 


