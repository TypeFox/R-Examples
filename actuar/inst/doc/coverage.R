### R code from vignette source 'coverage.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: coverage.Rnw:42-43
###################################################
library(actuar)


###################################################
### code chunk number 2: coverage.Rnw:88-90
###################################################
deductible <- 5
limit <- 13


###################################################
### code chunk number 3: coverage.Rnw:95-106
###################################################
pgammaL <- coverage(cdf = pgamma, deductible = deductible, limit = limit,
                     per.loss = TRUE)
dgammaL <- coverage(dgamma, pgamma, deductible = deductible, limit = limit,
                     per.loss = TRUE)
pgammaP <- coverage(cdf = pgamma, deductible = deductible, limit = limit)
dgammaP <- coverage(dgamma, pgamma, deductible = deductible, limit = limit)

d <- deductible
u <- limit - d
e <- 0.001
ylim <- c(0, dgammaL(0, 5, 0.6))


###################################################
### code chunk number 4: coverage.Rnw:131-137
###################################################
par(mar = c(2, 3, 1, 1))
curve(pgammaP(x, 5, 0.6), from = 0, to = u - e,
      xlim = c(0, limit), ylim = c(0, 1),
      xlab = "", ylab = "", xaxt = "n", lwd = 2)
curve(pgammaP(x, 5, 0.6), from = u, add = TRUE, lwd = 2)
axis(1, at = c(0, u), labels = c("0", "u - d"))


###################################################
### code chunk number 5: coverage.Rnw:154-160
###################################################
par(mar = c(2, 3, 1, 1))
curve(dgammaP(x, 5, 0.6), from = 0 + e, to = u - e,
      xlim = c(0, limit), ylim = ylim,
      xlab = "", ylab = "", xaxt = "n", lwd = 2)
points(u, dgammaP(u, 5, 0.6), pch = 16)
axis(1, at = c(0, u), labels = c("0", "u - d"))


###################################################
### code chunk number 6: coverage.Rnw:190-196
###################################################
par(mar = c(2, 3, 1, 1))
curve(pgammaL(x, 5, 0.6), from = 0, to = u - e,
      xlim = c(0, limit), ylim = c(0, 1),
      xlab = "", ylab = "", xaxt = "n", lwd = 2)
curve(pgammaL(x, 5, 0.6), from = u, add = TRUE, lwd = 2)
axis(1, at = c(0, u), labels = c("0", "u - d"))


###################################################
### code chunk number 7: coverage.Rnw:211-217
###################################################
par(mar = c(2, 3, 1, 1))
curve(dgammaL(x, 5, 0.6), from = 0 + e, to = u - e,
      xlim = c(0, limit), ylim = ylim,
      xlab = "", ylab = "", xaxt = "n", lwd = 2)
points(c(0, u), dgammaL(c(0, u), 5, 0.6), pch = 16)
axis(1, at = c(0, u), labels = c("0", "u - d"))


###################################################
### code chunk number 8: coverage.Rnw:225-238
###################################################
pgammaL <- coverage(cdf = pgamma, deductible = deductible, limit = limit,
                     per.loss = TRUE, franchise = TRUE)
dgammaL <- coverage(dgamma, pgamma, deductible = deductible, limit = limit,
                     per.loss = TRUE, franchise = TRUE)
pgammaP <- coverage(cdf = pgamma, deductible = deductible, limit = limit,
                    franchise = TRUE)
dgammaP <- coverage(dgamma, pgamma, deductible = deductible, limit = limit,
                    franchise = TRUE)

d <- deductible
u <- limit
e <- 0.001
ylim <- c(0, dgammaL(0, 5, 0.6))


###################################################
### code chunk number 9: coverage.Rnw:263-269
###################################################
par(mar = c(2, 3, 1, 1))
curve(pgammaP(x, 5, 0.6), from = 0, to = u - e,
      xlim = c(0, limit + d), ylim = c(0, 1),
      xlab = "", ylab = "", xaxt = "n", lwd = 2)
curve(pgammaP(x, 5, 0.6), from = u, add = TRUE, lwd = 2)
axis(1, at = c(0, d, u), labels = c("0", "d", "u"))


###################################################
### code chunk number 10: coverage.Rnw:286-293
###################################################
par(mar = c(2, 3, 1, 1))
curve(dgammaP(x, 5, 0.6), from = d + e, to = u - e,
      xlim = c(0, limit + d), ylim = ylim,
      xlab = "", ylab = "", xaxt = "n", lwd = 2)
curve(dgammaL(x, 5, 0.6), from = 0 + e, to = d, add = TRUE, lwd = 2)
points(u, dgammaP(u, 5, 0.6), pch = 16)
axis(1, at = c(0, d, u), labels = c("0", "d", "u"))


###################################################
### code chunk number 11: coverage.Rnw:323-329
###################################################
par(mar = c(2, 3, 1, 1))
curve(pgammaL(x, 5, 0.6), from = 0, to = u - e,
      xlim = c(0, limit + d), ylim = c(0, 1),
      xlab = "", ylab = "", xaxt = "n", lwd = 2)
curve(pgammaL(x, 5, 0.6), from = u, add = TRUE, lwd = 2)
axis(1, at = c(0, d, u), labels = c("0", "d", "u"))


###################################################
### code chunk number 12: coverage.Rnw:343-350
###################################################
par(mar = c(2, 3, 1, 1))
curve(dgammaL(x, 5, 0.6), from = d + e, to = u - e,
      xlim = c(0, limit + d), ylim = ylim,
      xlab = "", ylab = "", xaxt = "n", lwd = 2)
curve(dgammaL(x, 5, 0.6), from = 0 + e, to = d, add = TRUE, lwd = 2)
points(c(0, u), dgammaL(c(0, u), 5, 0.6), pch = 16)
axis(1, at = c(0, d, u), labels = c("0", "d", "u"))


