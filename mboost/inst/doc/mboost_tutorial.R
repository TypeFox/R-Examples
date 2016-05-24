### R code from vignette source 'mboost_tutorial.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: init
###################################################
## set up R session
options(prompt = "R> ", continue = "+  ", width = 80)
pd <- packageDescription("mboost")
if (any(is.na(pd))){
    install.packages("mboost", repos = "http://cran.at.r-project.org")
    pd <- packageDescription("mboost")
}
if (compareVersion(pd$Version, "2.1-0") < 0){ # must be mboost 2.1-X or newer
    warning("Current version of mboost is installed!")
    install.packages("mboost", repos = "http://cran.at.r-project.org")
}
require("mboost")
set.seed(190781)
### make graphics directory if not existing
if (!file.exists("graphics"))
    dir.create("graphics")


###################################################
### code chunk number 2: setup
###################################################
library("mboost")                         ## load package
data("bodyfat", package = "TH.data")      ## load data


###################################################
### code chunk number 3: mod_garcia
###################################################
## Reproduce formula of Garcia et al., 2005
lm1 <- lm(DEXfat ~ hipcirc  + kneebreadth + anthro3a, data = bodyfat)
coef(lm1)


###################################################
### code chunk number 4: mod_garcia_boosted
###################################################
## Estimate same model by glmboost
glm1 <- glmboost(DEXfat ~ hipcirc  + kneebreadth + anthro3a, data = bodyfat)
coef(glm1, off2int=TRUE)    ## off2int adds the offset to the intercept


###################################################
### code chunk number 5: mod_glmboost
###################################################
glm2 <- glmboost(DEXfat ~ ., data = bodyfat)


###################################################
### code chunk number 6: mod_glmboost2
###################################################
preds <- names(bodyfat[, names(bodyfat) != "DEXfat"])        ## names of predictors
fm <- as.formula(paste("DEXfat ~", paste(preds, collapse = "+")))  ## build formula
fm


###################################################
### code chunk number 7: coef_glmboost
###################################################
coef(glm2,        ## usually the argument 'which' is used to specify single base-
     which = "")  ## learners via partial matching; With which = "" we select all.


###################################################
### code chunk number 8: plot_glmboost
###################################################
plot(glm2, off2int = TRUE)    ## default plot, offset added to intercept
## now change ylim to the range of the coefficients without intercept (zoom-in)
plot(glm2, ylim = range(coef(glm2, which = preds)))


###################################################
### code chunk number 9: coefsplot1
###################################################
## same as first plot from chunk "plot_glmboost" but with better margins
par(mar = c(4, 4, 0, 6) + 0.1)
plot(glm2, main="", off2int=TRUE)


###################################################
### code chunk number 10: coefsplot2
###################################################
## same as second plot from chunk "plot_glmboost" but with better margins
par(mar = c(4, 4, 0, 6) + 0.1)
plot(glm2, ylim = range(coef(glm2, which = preds)), main = "")


###################################################
### code chunk number 11: mod_gamboost
###################################################
## now an additive model with the same variables as lm1
gam1 <- gamboost(DEXfat ~ bbs(hipcirc) + bbs(kneebreadth) + bbs(anthro3a),
                 data = bodyfat)


###################################################
### code chunk number 12: plot_gamboost
###################################################
par(mfrow = c(1,3))    ## 3 plots in one device
plot(gam1)             ## get the partial effects


###################################################
### code chunk number 13: partialplots
###################################################
## same as plot from chunk "plot_gamboost" but with better margins
par(mfrow=c(1,3))
par(mar = c(4, 4, 0, 0) + 0.1)
plot(gam1)


###################################################
### code chunk number 14: mod_gamboost2
###################################################
## every predictor enters the model via a bbs() base-learner (i.e., as smooth effect)
gam2 <- gamboost(DEXfat ~ ., baselearner = "bbs", data = bodyfat,
                 control = boost_control(trace = TRUE))


###################################################
### code chunk number 15: mod_gamboost2_refit
###################################################
## refit model to supress trace in cvrisk
## this is only required for the output in the PDF
gam2 <- gamboost(DEXfat ~ ., baselearner = "bbs", data = bodyfat)


###################################################
### code chunk number 16: mod_gamboost2_ctd
###################################################
set.seed(123)          ## set seed to make results reproducible
cvm <- cvrisk(gam2)    ## default method is 25-fold bootstrap cross-validation
## if package 'multicore' is not available this will trigger a warning
cvm
plot(cvm)    ## get the paths


###################################################
### code chunk number 17: cvpaths
###################################################
## same as plot from chunk "mod_gamboost2_ctd"with better margins
par(mar = c(4, 4, 0, 0) + 0.1)
plot(cvm, main = "")


###################################################
### code chunk number 18: mstop
###################################################
mstop(cvm)    ## extract the optimal mstop


###################################################
### code chunk number 19: mstop_ctd
###################################################
gam2[ mstop(cvm) ]    ## set the model automatically to the optimal mstop


###################################################
### code chunk number 20: mod_gamboost2_refit2
###################################################
## refit model to get trace again
## this is only required for the output in the PDF
gam2 <- gamboost(DEXfat ~ ., baselearner = "bbs", data = bodyfat,
                 control = boost_control(trace = TRUE))
gam2[ mstop(cvm) ]


###################################################
### code chunk number 21: coef_gamboost
###################################################
names(coef(gam2)) ## displays the selected base-learners at iteration 30
## To see that nothing got lost we now increase mstop to 1000:
gam2[1000, return = FALSE]  # return = FALSE just supresses "print(gam2)"


###################################################
### code chunk number 22: coef_gamboost2
###################################################
names(coef(gam2))  ## displays the selected base-learners, now at iteration 1000


###################################################
### code chunk number 23: quantreg_gamboost
###################################################
## Same model as glm1 but now with QuantReg() family
glm3 <- glmboost(DEXfat ~ hipcirc + kneebreadth + anthro3a, data = bodyfat,
                 family = QuantReg(tau = 0.5), control = boost_control(mstop = 500))
coef(glm3, off2int = TRUE)


###################################################
### code chunk number 24: loss
###################################################
loss =  function(y, f) tau * (y - f) * ((y - f) >= 0) +
                       (tau - 1) * (y - f) * ((y - f) < 0)


###################################################
### code chunk number 25: ngradient
###################################################
ngradient = function(y, f, w = NULL)  tau * ((y - f) >= 0) +
                    (tau- 1) * ((y - f) < 0)


###################################################
### code chunk number 26: OurQR
###################################################
OurQuantReg <- function(tau = 0.5){    ## function to include dependency on tau
  Family(                              ## applying the Family function
         loss = function(y, f)                               ## loss as above
                    tau       * (y - f) * ((y - f) >= 0) +
                    (tau - 1) * (y - f) * ((y - f) < 0) ,
         ngradient = function(y, f, w = NULL)                ## ngradient as above
                    tau * ((y - f) >= 0) + (tau - 1) * ((y - f) < 0),
         offset = function(y, w = NULL)                      ## median as offset
                    quantile(y, p = 0.5),
         name = "Our new family for quantile regression" )}
OurQuantReg()


###################################################
### code chunk number 27: ourQR_gamboost
###################################################
## Same model as glm3 but now with our new family
glm3b <- glmboost(DEXfat ~ hipcirc + kneebreadth + anthro3a, data = bodyfat,
                  family = OurQuantReg(tau = 0.5),
                  control = boost_control(mstop = 500))
identical(coef(glm3b), coef(glm3))


###################################################
### code chunk number 28: QR_hipcirc
###################################################
glm4a <- glmboost(DEXfat ~ hipcirc, family = OurQuantReg(tau = 0.05), data = bodyfat,
                  control = boost_control(mstop = 2000))
glm4b <- glmboost(DEXfat ~ hipcirc, family = OurQuantReg(tau = 0.5),  data = bodyfat,
                  control = boost_control(mstop = 2000))
glm4c <- glmboost(DEXfat ~ hipcirc, family = OurQuantReg(tau = 0.95), data = bodyfat,
                  control = boost_control(mstop = 2000))


###################################################
### code chunk number 29: plot_QR_hipcirc
###################################################
ord <- order(bodyfat$hipcirc)    ## order the data to avoid problems when plotting
plot(bodyfat$hipcirc[ord], bodyfat$DEXfat[ord])                   ## observed data
lines(bodyfat$hipcirc[ord], fitted(glm4a)[ord], lty = 2, lwd = 2) ## 0.05 quantile
lines(bodyfat$hipcirc[ord], fitted(glm4b)[ord], lty = 1, lwd = 2) ## median
lines(bodyfat$hipcirc[ord], fitted(glm4c)[ord], lty = 2, lwd = 2) ## 0.95 quantile


###################################################
### code chunk number 30: quantregbodyfat
###################################################
## same plot but with better margings and parameters
par(mar = c(4, 4, 0, 0) + 0.1)
plot(bodyfat$hipcirc[ord], bodyfat$DEXfat[ord])
lines(bodyfat$hipcirc[ord], fitted(glm4a)[ord], lty=2,lwd=2)
lines(bodyfat$hipcirc[ord], fitted(glm4b)[ord], lty=1, lwd=2)
lines(bodyfat$hipcirc[ord], fitted(glm4c)[ord], lty=2, lwd=2)


###################################################
### code chunk number 31: init_figs
###################################################
################################################################################
## the following chunks produce the graphics that are NOT part of the bodyfat ##
## example                                                                    ##
################################################################################


###################################################
### code chunk number 32: init
###################################################
red <- rgb(103,0,31, max = 255) ## define red color


###################################################
### code chunk number 33: center_false
###################################################
## load library mboost
library("mboost")
library("RColorBrewer")
## define colors
cols <- paste(brewer.pal(11,"RdBu")[c(10,2)], "E6", sep="")

## create graphics to show importance of centering
set.seed(1907)
x <- rnorm(50, mean = 5)
y <- rnorm(50, mean = x, sd = 0.3)
## subtract offset as usual in boosting
y <- y - mean(y)

par(ps = 8, cex=1, cex.lab=1, mar=c(3.1,3,0.5,0.1), mgp=c(2,1,0))
xrange <- range(0,x)
plot(y ~ x, xlim = xrange, cex = 1,
     xlab = "x", ylab = "negative gradient in step m = 1",
     pch = 20,
     col = rgb(0.5, 0.5, 0.5, alpha = 0.8))
abline(h = 0, col = "gray", lwd = 0.5)
abline(v = 0, col = "gray", lwd = 0.5)
abline(lm(y ~ x -1), col = cols[1], lwd = 1.5)
points(0, 0, col = cols[2], lwd = 1.5)
points(mean(x), mean(y), col = cols[2], pch = 3, lwd = 1.5)
legend(0.1, 2.35,  legend = c("origin", "center of data", "base-learner"),
       pch = c(21, 3, -1), lty = c(-1, -1, 1),
       col = c(cols[2], cols[2], cols[1]), lwd = 1.5, bg = "white", bty = "n")


###################################################
### code chunk number 34: center_true
###################################################
## create graphics to show importance of centering
set.seed(1907)
x <- rnorm(50, mean = 5)
y <- rnorm(50, mean = x, sd = 0.3)
## subtract offset as usual in boosting
y <- y - mean(y)

## figure with centering of x
mx <- mean(x)
xrange <- range(0,x)
x <- x - mean(x)

par(ps = 8, cex=1, cex.lab=1, mar=c(3.1,3,0.5,0.1), mgp=c(2,1,0))
plot(y ~ x, xlim = xrange - mx, cex = 1,
     xlab = "x (centered)", ylab = "negative gradient in step m = 1",
     pch = 20,
     col = rgb(0.5, 0.5, 0.5, alpha = 0.8))
abline(h = 0, col = "gray", lwd = 0.5)
abline(v = 0, col = "gray", lwd = 0.5)
abline(lm(y ~ x -1), col = cols[1], lwd = 1.5)
points(0, 0, col = cols[2], lwd = 1.5)
points(mean(x), mean(y), col = cols[2], pch = 3, lwd = 1.5)


###################################################
### code chunk number 35: bolsx1
###################################################
### Simulate some data
set.seed(1907)
n <- 100
x1 <- rnorm(n)
x2 <- as.factor(sample(0:3, n, replace = TRUE))
y <- 0.5 * x1 +  rnorm(n)
mod <- gamboost(y ~ bols(x1), control = boost_control(mstop = 25))
par(mar = c(4, 4, 0, 0) + 0.1)
plot(sort(x1), (0.5 * x1)[order(x1)], type = "l", lwd = 2,
     xlab = expression(x[1]), ylab = expression(f(x[1])))
lines(sort(x1), fitted(mod, which = 1)[order(x1)], col = red,
      lwd = 2, type = "b", pch = 20)
legend("topleft", c("true effect", "model"),
       lty = c(1, 1), pch = c(-1, 20), merge = TRUE,
       lwd = 2,
       col = c("black", red), bty = "n")


###################################################
### code chunk number 36: bolsx2
###################################################
beta <- c(0, -1, 0.5, 3)
y <- drop(model.matrix(~ x2) %*% beta + rnorm(n, sd = 0.3))
mod <- gamboost(y ~ bols(x2), control = boost_control(mstop = 50))
par(mar = c(4, 4, 0, 0) + 0.1)
betaPred <- coef(mod)[[1]]
betaPred[1] <- betaPred[1] + attr(coef(mod), "offset")
betaTrue <- c(0, -1, 0.5, 3)
plot(1:4, betaPred, type = "n", xaxt = "n",
     xlab = expression(x[2]), ylab = expression(f(x[2])),
     xlim = c(0.5, 4.5), ylim = c(-1.5, 3.5))
axis(1, at = 1:4, labels = expression(x[2]^(1), x[2]^(2), x[2]^(3), x[2]^(4)))
for (i in 1:4)
    lines(i + c(-0.38, 0, 0.38), rep(betaPred[i], each = 3),
          lwd = 2, col = red, type = "b", pch = 20)
for (i in 1:4)
    lines(i + c(-0.4, 0.4), rep(betaTrue[i], each = 2),
          lwd = 2, col = "black")
legend("topleft", c("true effect", "model"),
       pch = c(-1,20), lty = c(1, 1), lwd = 2,
       col = c("black", red), bty = "n")


###################################################
### code chunk number 37: bbsx1
###################################################
set.seed(1907)
n <- 100
x1 <- rnorm(n)
x2 <- rnorm(n) + 0.25 * x1
x3 <- as.factor(sample(0:1, n, replace = TRUE))
x4 <- gl(4, 25)
y <- 3 * sin(x1) + x2^2 + rnorm(n)

## specify knot locations for bbs(x2, ...)
knots.x2 <- quantile(x2, c(0.25, 0.5, 0.75))

## fit the model
mod <- gamboost(y ~ bbs(x1, knots = 20, df = 4) +
                bbs(x2, knots = knots.x2, df = 5) +
                bols(x3) + bols(x4))

## plot partial effects
par(mar = c(4, 4, 0, 0) + 0.1)
plot(sort(x1), (3 * sin(x1))[order(x1)], type = "l", lwd = 2,
     xlab = expression(x[1]), ylab =  expression(f[1](x[1])))
lines(sort(x1), fitted(mod, which = 1)[order(x1)], col = red, lwd = 2,
      type = "b", pch = 20)
legend("topleft", c("true effect", "model"),
       lty = c(1, 1), pch = c(-1,20), merge = TRUE, lwd = 2,
       col = c("black", red), bty = "n")


###################################################
### code chunk number 38: bbsx2
###################################################
par(mar = c(4, 4, 0, 0) + 0.1)
plot(sort(x2), (x2^2)[order(x2)], type = "l", lwd = 2,
     xlab = expression(x[2]), ylab =  expression(f[2](x[2])))
## offset needed to conserve correct level (otherwise center both effects around
## zero to make them comparable)
lines(sort(x2), fitted(mod, which = 2)[order(x2)] + mod$offset, col = red, lwd = 2,
      type = "b", pch = 20)


###################################################
### code chunk number 39: cyclic1
###################################################
### example for cyclic spline
set.seed(1907)
true_f <- function(x)
    cos(x) + 0.25 * sin(4 * x)
x <- runif(150, 0, 2 * pi)
y <- rnorm(150, mean = true_f(x), sd=0.1)
newX <- seq(0, 2*pi, length = 200)
mod <- gamboost(y ~ bbs(x, knots = 20, degree = 4))
mod[3000]
mod_cyclic <- gamboost(y ~ bbs(x, cyclic=TRUE, knots = 20,
                               degree = 4, boundary.knots=c(0, 2*pi)))
mod_cyclic[3000]

par(mar = c(4, 4, 0, 0) + 0.1)
plot(x,y, ylab = "f(x)",
     pch = 20, xlim = c(0, 7),
     col = rgb(0.5, 0.5, 0.5, alpha = 0.8))
lines(newX, predict(mod, data.frame(x = newX)),
      col = "black", lwd = 2)
lines(newX + 2 * pi, predict(mod, data.frame(x = newX)),
      col = "black", lwd = 2)
legend("bottomleft", c("cyclic = FALSE"),
       lty = c(1), col = c("black"), lwd = 2, bty = "n")


###################################################
### code chunk number 40: cyclic2
###################################################
par(mar = c(4, 4, 0, 0) + 0.1)
plot(x, y, ylab = "f(x)",
     pch = 20, xlim = c(0, 7),
     col = rgb(0.5, 0.5, 0.5, alpha = 0.8))
lines(newX, predict(mod_cyclic, data.frame(x = newX)),
      col = red, lwd = 2)
lines(newX + 2 * pi, predict(mod_cyclic, data.frame(x = newX)),
      col = red, lwd = 2)
abline(v = 2 * pi, col = "gray", lty = "dashed", lwd = 2)
legend("bottomleft", c("cyclic = TRUE"),
       lty = c(1), col = c(red), lwd = 2, bty = "n")


###################################################
### code chunk number 41: bspatial1
###################################################
set.seed(1907)
x1 <- runif(250,-pi,pi)
x2 <- runif(250,-pi,pi)
y <- sin(x1) * sin(x2) + rnorm(250, sd = 0.4)
modspa <- gamboost(y ~ bspatial(x1, x2, knots = list(x1 = 12, x2 = 12)))
# possible to specify different knot mesh for x1 and x2:
# modspa <- gamboost(y ~ bspatial(x1, x2, knots = list(x1 = 12, x2 = 20)))

library(lattice)
library(RColorBrewer)
## set up colorscheme
nCols <- 50
cols <- colorRampPalette(rev(brewer.pal(11, "RdBu")), space = "Lab",
                         interpolate = "spline")(nCols)
## make new data on a grid:
nd <- expand.grid(x1 = seq(-pi, pi, length = 100),
                  x2 = seq(-pi, pi, length = 100))
preds <- predict(modspa, newdata = nd)
breaks <- seq(min(preds), max(preds), length = nCols + 1)
print(levelplot(preds ~ nd$x1 + nd$x2,
          xlab = expression(x[1]),
          ylab = expression(x[2]),
          at = breaks,
          col.regions = cols, cex = 0.7,
          colorkey = list(space = "left", at = breaks, width = 1)))


###################################################
### code chunk number 42: bspatial2
###################################################
x1 <- x2 <- seq(-pi, pi, length = 50)
nd <- expand.grid(x1 = x1,
                  x2 = x2)
preds <- predict(modspa, newdata = nd)
z <- matrix(preds, ncol = length(x1))
nrz <- nrow(z)
ncz <- ncol(z)
jet.colors <- colorRampPalette(paste(brewer.pal(11,"RdBu")[c(10,2)],
                                     "E6", sep=""))
nbcol <- 10
color <- jet.colors(nbcol)
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)
par(mar=c(1, 0.1, 0.1, 0.1), mgp = c(1.8,1,0))
persp(x1, x2, z,
      theta = 45, phi = 30, expand = 0.5,
      col = color[facetcol], ticktype = "detailed",
      nticks = 3,
      xlab = "x1", ylab = "x2")


###################################################
### code chunk number 43: families
###################################################
pdf("./graphics/fig-family.pdf", width = 5, height = 4)
par(mar = c(4, 4, 0, 0) + 0.1)
## Gaussian
x <- seq(-2, 2, length = 200)
plot(x,x^2, type = "l", ylab = expression(rho), xlab = expression(y-f))

## Laplace
x <- seq(-2, 2, length = 200)
lossL1 <- function(x, param) abs(x)
mp <- 0.5
dat <- data.frame(x = rep(x, length(mp)),
                  y = as.numeric(sapply(mp, function(z) lossL1(x, param = z))),
                  param = rep(mp, each = length(x)))
dat$tau <- factor(dat$param)
plot(x, dat$y, type = "l", ylab = expression(rho), xlab = expression(y-f))

## Huber
x <- seq(-2, 2, length = 200)
lossH <- function(x, param) ifelse(abs(x) <= param, (1/2) * x^2, param * (abs(x) - param / 2))
mp <- c(0.2, 0.5, 1, 2, 10)
dat <- data.frame(x = x, sapply(mp, function(z) lossH(x, param = z)))
plot(x, dat$X1, type = "l", ylim = c(0,2), ylab = expression(rho),
     xlab = expression(y-f))
legend(0, 2, xjust = 0.5, legend = paste(mp), title = expression(delta),
       lty = 1, col = c(1,3,4,5,6), cex = 0.6, box.col = "gray")
for(i in 1:(length(mp)-1)){
    lines(x, dat[,i+2], col = i+2)
}

## Quantile
x <- seq(-2, 2, length = 200)
lossL1 <- function(x, param) ifelse(x >= 0, param  * x, (param - 1) * x)
mp <- c(5:9 / 10)
dat <- data.frame(x = x, sapply(mp, function(z) lossL1(x, param = z)))
plot(x, dat$X1, type = "l", ylab = expression(rho), ylim = c(0,2),
     xlab = expression(y-f))
legend(0,2, xjust = 0.5, legend = paste(mp), title = expression(tau),
       lty = 1, col = c(1,3,4,5,6), cex = 0.6, box.col = "gray")
for(i in 1:(length(mp)-1)){
    lines(x, dat[,i+2], col = i+2)
}

## Binary
x <- seq(-2, 2, length = 200)
dat <- data.frame(x = x, Binomial = log(1+exp(-2*x), base=2), AdaExp = exp(-x))
plot(x, dat$Binomial, type = "l", lty = 1, col = 1, xlab = expression(tilde(y) * f),
     ylab = expression(rho), ylim = c(0,7.5))
lines(x, dat$AdaExp, lty = 2, col = 2)
legend("topright", legend = c("Binomial","AdaExp"), xjust = 0.5, title = "loss",
       lty = 1:2, col = 1:2, cex = 0.8, bty = "n")
dev.off()


