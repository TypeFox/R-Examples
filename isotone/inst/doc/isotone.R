### R code from vignette source 'isotone.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("isotone")
options(prompt = "R> ", continue = "+  ")


###################################################
### code chunk number 2: isotone.Rnw:447-450
###################################################
require("isotone")
data("pituitary")
head(pituitary)


###################################################
### code chunk number 3: isotone.Rnw:455-458
###################################################
res1 <- with(pituitary, gpava(age, size, ties = "primary"))
res2 <- with(pituitary, gpava(age, size, ties = "secondary"))
res3 <- with(pituitary, gpava(age, size, ties = "tertiary"))


###################################################
### code chunk number 4: pit-plot (eval = FALSE)
###################################################
## layout(matrix(c(1,1,2,2,0,3,3,0), 2, 4, byrow = TRUE))
## plot(res1, main = "PAVA plot (primary)")
## plot(res2, main = "PAVA plot (secondary)")
## plot(res3,  main = "PAVA plot (tertiary)")


###################################################
### code chunk number 5: pit-plot1
###################################################
layout(matrix(c(1,1,2,2,0,3,3,0), 2, 4, byrow = TRUE))
plot(res1, main = "PAVA plot (primary)")
plot(res2, main = "PAVA plot (secondary)")
plot(res3,  main = "PAVA plot (tertiary)")


###################################################
### code chunk number 6: isotone.Rnw:480-481
###################################################
tapply(res3$x, res3$z, mean)


###################################################
### code chunk number 7: isotone.Rnw:493-495
###################################################
data("posturo")
head(posturo)


###################################################
### code chunk number 8: isotone.Rnw:500-504
###################################################
res.mean <- with(posturo, gpava(height, cbind(SOT.1, SOT.2, SOT.3),
  solver = weighted.mean, ties = "secondary"))
res.median <- with(posturo, gpava(height, cbind(SOT.1, SOT.2, SOT.3),
  solver = weighted.median, ties = "secondary"))


###################################################
### code chunk number 9: post-plot (eval = FALSE)
###################################################
## plot(res.mean)
## plot(res.median)


###################################################
### code chunk number 10: post-plot1
###################################################
par(mar = c(5,4,4,2))
par(mfrow = c(1, 2))
plot(res.mean)
plot(res.median)


###################################################
### code chunk number 11: isotone.Rnw:540-544
###################################################
set.seed(12345)
y <- rnorm(9)
w1 <- rep(1, 9)
Atot <- cbind(1:8, 2:9)


###################################################
### code chunk number 12: isotone.Rnw:549-550
###################################################
Atot


###################################################
### code chunk number 13: isotone.Rnw:557-561
###################################################
fit.ls1 <- activeSet(Atot, "LS", y = y, weights = w1)
fit.ls2 <- activeSet(Atot, fSolver, y = y, weights = w1,
  fobj = function(x) sum(w1 * (x - y)^2),
  gobj = function(x) 2 * drop(w1 * (x - y)))


###################################################
### code chunk number 14: isotone.Rnw:566-571
###################################################
set.seed(12345)
wvec <- 1:9
wmat <- crossprod(matrix(rnorm(81), 9, 9))/9
fit.wls <- activeSet(Atot, "LS", y = y, weights = wvec)
fit.gls <- activeSet(Atot, "GLS", y = y, weights = wmat)


###################################################
### code chunk number 15: isotone.Rnw:575-576
###################################################
fit.qua <- activeSet(Atot, "quantile", y = y, weights = wvec, aw = 0.3, bw = 0.7)


###################################################
### code chunk number 16: isotone.Rnw:581-582
###################################################
fit.abs <- activeSet(Atot, "L1", y = y, weights = w1)


###################################################
### code chunk number 17: isotone.Rnw:587-589
###################################################
fit.eps <- activeSet(Atot, "L1eps", y = y, weights = w1, eps = 1e-04)
fit.pow <- activeSet(Atot, "Lp", y = y, weights = w1, p = 1.2)


###################################################
### code chunk number 18: isotone.Rnw:594-595
###################################################
fit.che <- activeSet(Atot, "chebyshev", y = y, weights = w1)


###################################################
### code chunk number 19: isotone.Rnw:600-601
###################################################
fit.asy <- activeSet(Atot, "asyLS", y = y, weights = w1, aw = 2, bw = 1)


###################################################
### code chunk number 20: isotone.Rnw:606-608
###################################################
fit.hub <- activeSet(Atot, "huber", y = y, weights = w1, eps = 1)
fit.svm <- activeSet(Atot, "SILF", y = y, weights = w1, beta = 0.8, eps = 0.2)


###################################################
### code chunk number 21: isotone.Rnw:613-617
###################################################
set.seed(12345)
yp <- rpois(9, 5)
x0 <- 1:9
fit.poi <- activeSet(Atot, "poisson", x0 = x0, y = yp)


###################################################
### code chunk number 22: isotone.Rnw:622-625
###################################################
Atree <- matrix(c(1, 1, 2, 2, 2, 3, 3, 8, 2, 3, 4, 5, 6, 7, 8, 9), 8, 2)
Atree
fit.tree <- activeSet(Atree, "LS", y = y, weights = w1)


###################################################
### code chunk number 23: isotone.Rnw:630-634
###################################################
Aloop <- matrix(c(1, 2, 3, 3, 4, 5, 6, 6, 7, 8, 3, 3, 4, 5, 6, 6,
  7, 8, 9, 9), 10, 2)
Aloop
fit.loop <- activeSet(Aloop, "LS", y = y, weights = w1)


###################################################
### code chunk number 24: isotone.Rnw:639-643
###################################################
Ablock <- cbind(c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3),
  rep(6, 3)), c(rep(c(4, 5, 6), 3), rep(c(7, 8, 9), 3)))
Ablock
fit.block <- activeSet(Ablock, "LS", y = y, weights = w1)


###################################################
### code chunk number 25: isotone.Rnw:651-655
###################################################
pava.fitted <- gpava(1:9, y)$x
aset.fitted <- activeSet(Atot, "LS", weights = w1, y = y)$x
mse <- mean((pava.fitted - aset.fitted)^2)
mse


