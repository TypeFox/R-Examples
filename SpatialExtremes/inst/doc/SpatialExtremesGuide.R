### R code from vignette source 'SpatialExtremesGuide.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: SpatialExtremesGuide.Rnw:32-42
###################################################
getOption("SweaveHooks")[["fig"]]()
library(SpatialExtremes)
set.seed(12)
n.site <- 50
coord <- matrix(runif(2*n.site, 0, 10), ncol = 2)

##Simulate a max-stable process - with unit Frechet margins
data <- rmaxstab(40, coord, "whitmat", nugget = 0, range = 1, smooth = 2)

##Compute the lambda-madogram
lmadogram(data, coord, n.bins = 25, border = "grey", box = FALSE)


###################################################
### code chunk number 2: SpatialExtremesGuide.Rnw:55-57
###################################################
options(SweaveHooks=list(fig=function()
par(mar=c(5.1, 4.1, 1.1, 2.1))))


###################################################
### code chunk number 3: SpatialExtremesGuide.Rnw:98-99 (eval = FALSE)
###################################################
## vignette("SpatialExtremesGuide")


###################################################
### code chunk number 4: Smith2Sim (eval = FALSE)
###################################################
## x <- seq(0, 10, length = 100)
## seed <- 19
## set.seed(seed)
## data0 <- rmaxstab(1, cbind(x,x), "gauss", cov11 = 9/8, cov12 = 0, cov22 = 9/8,
##                   grid = TRUE)
## set.seed(seed)
## data1 <- rmaxstab(1, cbind(x,x), "gauss", cov11 = 9/8, cov12 = 1, cov22 = 9/8,
##                   grid = TRUE)
## 
## png("Figures/Smith2Sim.png", width = 1400, height = 700)
## par(mfrow=c(1,2))
## image(x, x, log(data0), col = terrain.colors(64))
## image(x, x, log(data1), col = terrain.colors(64))
## dev.off()


###################################################
### code chunk number 5: covariances
###################################################
par(mfrow=c(1,4), mar = c(4, 4, 0.2, 0.2))
covariance(nugget = 0, sill = 1, range = 1, smooth = 4, cov.mod = "whitmat",
           xlim = c(0,9), ylim = c(0, 1))
covariance(nugget = 0, sill = 1, range = 1, smooth = 2, cov.mod = "whitmat",
           add = TRUE, col = 2, xlim = c(0,9))
covariance(nugget = 0, sill = 1, range = 1, smooth = 1, cov.mod = "whitmat",
           add = TRUE, col = 3, xlim = c(0,9))
covariance(nugget = 0, sill = 1, range = 1, smooth = 0.5, cov.mod = "whitmat",
           col = 4, add = TRUE, xlim = c(0,9))
legend("topright", c(expression(nu == 4), expression(nu == 2),
                     expression(nu == 1), expression(nu == 0.5)),
       col = 1:4, lty = 1, inset = 0.05)

covariance(nugget = 0, sill = 1, range = 1, smooth = 2, cov.mod = "powexp",
           xlim = c(0, 4), ylim = c(0, 1))
covariance(nugget = 0, sill = 1, range = 1, smooth = 1.5, cov.mod = "powexp",
           add = TRUE, col = 2, xlim = c(0, 4))
covariance(nugget = 0, sill = 1, range = 1, smooth = 1, cov.mod = "powexp",
           add = TRUE, col = 3, xlim = c(0, 4))
covariance(nugget = 0, sill = 1, range = 1, smooth = 0.75, cov.mod = "powexp",
           add = TRUE, col = 4, xlim = c(0, 4))
legend("topright", c(expression(nu == 2), expression(nu == 1.5),
                     expression(nu == 1), expression(nu == 0.75)),
       col = 1:4, lty = 1, inset = 0.05)

covariance(nugget = 0, sill = 1, range = 1, smooth = 4, cov.mod = "cauchy",
           xlim = c(0, 5), ylim = c(0, 1))
covariance(nugget = 0, sill = 1, range = 1, smooth = 2, cov.mod = "cauchy",
           add = TRUE, col = 2, xlim = c(0, 5))
covariance(nugget = 0, sill = 1, range = 1, smooth = 1, cov.mod = "cauchy",
           add = TRUE, col = 3, xlim = c(0, 5))
covariance(nugget = 0, sill = 1, range = 1, smooth = 0.75, cov.mod = "cauchy",
           add = TRUE, col = 4, xlim = c(0, 5))
legend("topright", c(expression(nu == 4), expression(nu == 2),
                     expression(nu == 1), expression(nu == 0.75)),
       col = 1:4, lty = 1, inset = 0.05)

covariance(nugget = 0, sill = 1, range = 1, smooth = 2, cov.mod = "bessel",
           xlim = c(0, 20), ylim = c(-0.3, 1))
covariance(nugget = 0, sill = 1, range = 1, smooth = 1, cov.mod = "bessel",
           add = TRUE, col = 2, xlim = c(0, 20))
covariance(nugget = 0, sill = 1, range = 1, smooth = .5, cov.mod = "bessel",
           add = TRUE, col = 3, xlim = c(0, 20))
covariance(nugget = 0, sill = 1, range = 1, smooth = 0.25, cov.mod = "bessel",
           add = TRUE, col = 4, xlim = c(0, 20))
legend("topright", c(expression(nu == 2), expression(nu == 1),
                     expression(nu == .5), expression(nu == 0.25)),
       col = 1:4, lty = 1, inset = 0.05)


###################################################
### code chunk number 6: SpatialExtremesGuide.Rnw:554-555
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(1,4), mar = c(4, 4, 0.2, 0.2))
covariance(nugget = 0, sill = 1, range = 1, smooth = 4, cov.mod = "whitmat",
           xlim = c(0,9), ylim = c(0, 1))
covariance(nugget = 0, sill = 1, range = 1, smooth = 2, cov.mod = "whitmat",
           add = TRUE, col = 2, xlim = c(0,9))
covariance(nugget = 0, sill = 1, range = 1, smooth = 1, cov.mod = "whitmat",
           add = TRUE, col = 3, xlim = c(0,9))
covariance(nugget = 0, sill = 1, range = 1, smooth = 0.5, cov.mod = "whitmat",
           col = 4, add = TRUE, xlim = c(0,9))
legend("topright", c(expression(nu == 4), expression(nu == 2),
                     expression(nu == 1), expression(nu == 0.5)),
       col = 1:4, lty = 1, inset = 0.05)

covariance(nugget = 0, sill = 1, range = 1, smooth = 2, cov.mod = "powexp",
           xlim = c(0, 4), ylim = c(0, 1))
covariance(nugget = 0, sill = 1, range = 1, smooth = 1.5, cov.mod = "powexp",
           add = TRUE, col = 2, xlim = c(0, 4))
covariance(nugget = 0, sill = 1, range = 1, smooth = 1, cov.mod = "powexp",
           add = TRUE, col = 3, xlim = c(0, 4))
covariance(nugget = 0, sill = 1, range = 1, smooth = 0.75, cov.mod = "powexp",
           add = TRUE, col = 4, xlim = c(0, 4))
legend("topright", c(expression(nu == 2), expression(nu == 1.5),
                     expression(nu == 1), expression(nu == 0.75)),
       col = 1:4, lty = 1, inset = 0.05)

covariance(nugget = 0, sill = 1, range = 1, smooth = 4, cov.mod = "cauchy",
           xlim = c(0, 5), ylim = c(0, 1))
covariance(nugget = 0, sill = 1, range = 1, smooth = 2, cov.mod = "cauchy",
           add = TRUE, col = 2, xlim = c(0, 5))
covariance(nugget = 0, sill = 1, range = 1, smooth = 1, cov.mod = "cauchy",
           add = TRUE, col = 3, xlim = c(0, 5))
covariance(nugget = 0, sill = 1, range = 1, smooth = 0.75, cov.mod = "cauchy",
           add = TRUE, col = 4, xlim = c(0, 5))
legend("topright", c(expression(nu == 4), expression(nu == 2),
                     expression(nu == 1), expression(nu == 0.75)),
       col = 1:4, lty = 1, inset = 0.05)

covariance(nugget = 0, sill = 1, range = 1, smooth = 2, cov.mod = "bessel",
           xlim = c(0, 20), ylim = c(-0.3, 1))
covariance(nugget = 0, sill = 1, range = 1, smooth = 1, cov.mod = "bessel",
           add = TRUE, col = 2, xlim = c(0, 20))
covariance(nugget = 0, sill = 1, range = 1, smooth = .5, cov.mod = "bessel",
           add = TRUE, col = 3, xlim = c(0, 20))
covariance(nugget = 0, sill = 1, range = 1, smooth = 0.25, cov.mod = "bessel",
           add = TRUE, col = 4, xlim = c(0, 20))
legend("topright", c(expression(nu == 2), expression(nu == 1),
                     expression(nu == .5), expression(nu == 0.25)),
       col = 1:4, lty = 1, inset = 0.05)


###################################################
### code chunk number 7: SpatialExtremesGuide.Rnw:608-619 (eval = FALSE)
###################################################
## covariance(nugget = 0, sill = 1, range = 1, smooth = 4, cov.mod = "whitmat",
##            xlim = c(0,9), ylim = c(0, 1))
## covariance(nugget = 0, sill = 1, range = 1, smooth = 2, cov.mod = "whitmat",
##            add = TRUE, col = 2, xlim = c(0,9))
## covariance(nugget = 0, sill = 1, range = 1, smooth = 1, cov.mod = "whitmat",
##            add = TRUE, col = 3, xlim = c(0,9))
## covariance(nugget = 0, sill = 1, range = 1, smooth = 0.5, cov.mod = "whitmat",
##            col = 4, add = TRUE, xlim = c(0,9))
## legend("topright", c(expression(nu == 4), expression(nu == 2),
##                      expression(nu == 1), expression(nu == 0.5)),
##        col = 1:4, lty = 1, inset = 0.05)


###################################################
### code chunk number 8: Schlather2Sim (eval = FALSE)
###################################################
## x <- y <- seq(0, 10, length = 100)
## set.seed(12)
## ms0 <- rmaxstab(1, cbind(x, y), "whitmat", nugget = 0, range = 1, smooth = 1, grid=TRUE)
## set.seed(12)
## ms1 <- rmaxstab(1, cbind(x, y), "powexp", nugget = 0, range = 1.5, smooth = 1, grid=TRUE)
## 
## png("Figures/Schlather2Sim.png", width = 1400, height = 700)
## par(mfrow=c(1,2))
## image(x, y, log(ms0), col = terrain.colors(64))
## image(x, y, log(ms1), col = terrain.colors(64))
## dev.off()


###################################################
### code chunk number 9: anisoCovFun
###################################################
cov.fun <- covariance(nugget = 0, sill = 1, range = 1, smooth = 1,
                      cov.mod = "powexp", plot = FALSE)
phi <- pi / 4
r <- sqrt(0.5)
A <- matrix(c(cos(phi), r^2 * sin(phi), r^2 * sin(phi), cos(phi)), 2)
cov.fun.aniso <- function(vec.pos)
  cov.fun(sqrt(vec.pos %*% A %*% vec.pos))

rho1 <- rho2 <- matrix(NA, 100, 100)
xs <- ys <- seq(-4, 4, length = 100)
for (i in 1:100){
  x <- xs[i]
  for (j in 1:100){
    rho1[i,j] <- cov.fun(sqrt(x^2+ys[j]^2))
    rho2[i,j] <- cov.fun.aniso(c(x, ys[j]))
  }
}

par(mfrow=c(1,2))
contour(xs, ys, rho1, xlab = expression(paste(Delta, x)),
        ylab = expression(paste(Delta, y)))
contour(xs, ys, rho2, xlab = expression(paste(Delta, x)),
        ylab = expression(paste(Delta, y)))


###################################################
### code chunk number 10: SpatialExtremesGuide.Rnw:684-685
###################################################
getOption("SweaveHooks")[["fig"]]()
cov.fun <- covariance(nugget = 0, sill = 1, range = 1, smooth = 1,
                      cov.mod = "powexp", plot = FALSE)
phi <- pi / 4
r <- sqrt(0.5)
A <- matrix(c(cos(phi), r^2 * sin(phi), r^2 * sin(phi), cos(phi)), 2)
cov.fun.aniso <- function(vec.pos)
  cov.fun(sqrt(vec.pos %*% A %*% vec.pos))

rho1 <- rho2 <- matrix(NA, 100, 100)
xs <- ys <- seq(-4, 4, length = 100)
for (i in 1:100){
  x <- xs[i]
  for (j in 1:100){
    rho1[i,j] <- cov.fun(sqrt(x^2+ys[j]^2))
    rho2[i,j] <- cov.fun.aniso(c(x, ys[j]))
  }
}

par(mfrow=c(1,2))
contour(xs, ys, rho1, xlab = expression(paste(Delta, x)),
        ylab = expression(paste(Delta, y)))
contour(xs, ys, rho2, xlab = expression(paste(Delta, x)),
        ylab = expression(paste(Delta, y)))


###################################################
### code chunk number 11: tbm (eval = FALSE)
###################################################
## x <- y <- seq(0, 10, length = 100)
## coord <- cbind(x, y)
## seed <- 3
## set.seed(seed)
## data1 <- rgp(1, coord, cov.mod = "whitmat", sill = 1, range = 1, smooth = 1,
##              grid = TRUE, control = list(nlines = 1))
## set.seed(seed)
## data2 <- rgp(1, coord, cov.mod = "whitmat", sill = 1, range = 1, smooth = 1,
##              grid = TRUE, control = list(nlines = 25))
## set.seed(seed)
## data3 <- rgp(1, coord, cov.mod = "whitmat", sill = 1, range = 1, smooth = 1,
##              grid = TRUE, control = list(nlines = 625))
## 
## png("Figures/tbm.png", width = 1400, height = 700)
## par(mfrow=c(1,3))
## image(x, y, data1, col = terrain.colors(64))
## image(x, y, data2, col = terrain.colors(64))
## image(x, y, data3, col = terrain.colors(64))
## dev.off()


###################################################
### code chunk number 12: SpatialExtremesGuide.Rnw:1007-1008 (eval = FALSE)
###################################################
## x <- y <- seq(0, 10, length = 100)
## coord <- cbind(x, y)
## seed <- 3
## set.seed(seed)
## data1 <- rgp(1, coord, cov.mod = "whitmat", sill = 1, range = 1, smooth = 1,
##              grid = TRUE, control = list(nlines = 1))
## set.seed(seed)
## data2 <- rgp(1, coord, cov.mod = "whitmat", sill = 1, range = 1, smooth = 1,
##              grid = TRUE, control = list(nlines = 25))
## set.seed(seed)
## data3 <- rgp(1, coord, cov.mod = "whitmat", sill = 1, range = 1, smooth = 1,
##              grid = TRUE, control = list(nlines = 625))
## 
## png("Figures/tbm.png", width = 1400, height = 700)
## par(mfrow=c(1,3))
## image(x, y, data1, col = terrain.colors(64))
## image(x, y, data2, col = terrain.colors(64))
## image(x, y, data3, col = terrain.colors(64))
## dev.off()


###################################################
### code chunk number 13: rmaxstabSmith (eval = FALSE)
###################################################
## x <- y <- seq(0, 10, length = 100)
## coord <- cbind(x, y)
## set.seed(8)
## M0 <- rmaxstab(1, coord, "gauss", cov11 = 9/8, cov12 = 0, cov22 = 9/8,
##                grid = TRUE)
## set.seed(8)
## M1 <- rmaxstab(1, coord, "gauss", cov11 = 9/8, cov12 = 3/4, cov22 = 9/8,
##                grid = TRUE)
## 
## png("Figures/rmaxstabSmith.png", width = 1400, height = 700)
## par(mfrow = c(1,2))
## image(x, y, log(M0), col = terrain.colors(64))
## image(x, y, log(M1), col = terrain.colors(64))
## dev.off()


###################################################
### code chunk number 14: SpatialExtremesGuide.Rnw:1162-1163 (eval = FALSE)
###################################################
## x <- y <- seq(0, 10, length = 100)
## coord <- cbind(x, y)
## set.seed(8)
## M0 <- rmaxstab(1, coord, "gauss", cov11 = 9/8, cov12 = 0, cov22 = 9/8,
##                grid = TRUE)
## set.seed(8)
## M1 <- rmaxstab(1, coord, "gauss", cov11 = 9/8, cov12 = 3/4, cov22 = 9/8,
##                grid = TRUE)
## 
## png("Figures/rmaxstabSmith.png", width = 1400, height = 700)
## par(mfrow = c(1,2))
## image(x, y, log(M0), col = terrain.colors(64))
## image(x, y, log(M1), col = terrain.colors(64))
## dev.off()


###################################################
### code chunk number 15: SpatialExtremesGuide.Rnw:1189-1191 (eval = FALSE)
###################################################
## coord <- matrix(runif(100, 0, 10), ncol = 2)
## data1 <- rmaxstab(100, coord, "whitmat", nugget = 0, range = 1, smooth = 1)


###################################################
### code chunk number 16: SpatialExtremesGuide.Rnw:1194-1198 (eval = FALSE)
###################################################
## x <- seq(0, 10, length = 100)
## coord <- cbind(x, x)
## data2 <- rmaxstab(1, coord, "powexp", nugget = 0, range = 1, smooth = 2,
##                   grid = TRUE)


###################################################
### code chunk number 17: extCoeffModels
###################################################
smith <- function(h) 2 * pnorm(h/2)
cov.fun <- covariance(nugget = 0, sill = 1, range = 1, smooth = 1, plot = FALSE)
schlather <- function(h)
  1 + sqrt((1-cov.fun(h))/2)
alpha <- 0.1

plot(smith, from = 0, to = 5, xlab = "h", ylab = expression(theta(h)))
plot(schlather, add = TRUE, col = 2, from = 0, to = 5)
legend("bottomright", c("Smith", "Schlather"), col = 1:2, lty = 1, inset = 0.05)


###################################################
### code chunk number 18: SpatialExtremesGuide.Rnw:1303-1304
###################################################
getOption("SweaveHooks")[["fig"]]()
smith <- function(h) 2 * pnorm(h/2)
cov.fun <- covariance(nugget = 0, sill = 1, range = 1, smooth = 1, plot = FALSE)
schlather <- function(h)
  1 + sqrt((1-cov.fun(h))/2)
alpha <- 0.1

plot(smith, from = 0, to = 5, xlab = "h", ylab = expression(theta(h)))
plot(schlather, add = TRUE, col = 2, from = 0, to = 5)
legend("bottomright", c("Smith", "Schlather"), col = 1:2, lty = 1, inset = 0.05)


###################################################
### code chunk number 19: extCoeffST-Smith
###################################################
n.site <- 40
n.obs <- 100
coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
data <- rmaxstab(n.obs, coord, "whitmat", nugget = 0, range = 1,
                 smooth = 1)
cov.fun <- covariance(nugget = 0, sill = 1, range = 1, smooth = 1, cov.mod = "whitmat",
                      plot = FALSE)
ext.coeff <- function(h)
  1 + sqrt((1 - cov.fun(h))/2)

par(mfrow=c(1,2))
ST <- fitextcoeff(data, coord, loess = FALSE, ylim = c(1, 2))
plot(ext.coeff, from = 0, to = 12, col = 2, add = TRUE, lwd = 1.5)
Smith <- fitextcoeff(data, coord, estim = "Smith", loess = FALSE, ylim = c(1, 2))
plot(ext.coeff, from = 0, to = 12, col = 2, add = TRUE, lwd = 1.5)


###################################################
### code chunk number 20: SpatialExtremesGuide.Rnw:1407-1408
###################################################
getOption("SweaveHooks")[["fig"]]()
n.site <- 40
n.obs <- 100
coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
data <- rmaxstab(n.obs, coord, "whitmat", nugget = 0, range = 1,
                 smooth = 1)
cov.fun <- covariance(nugget = 0, sill = 1, range = 1, smooth = 1, cov.mod = "whitmat",
                      plot = FALSE)
ext.coeff <- function(h)
  1 + sqrt((1 - cov.fun(h))/2)

par(mfrow=c(1,2))
ST <- fitextcoeff(data, coord, loess = FALSE, ylim = c(1, 2))
plot(ext.coeff, from = 0, to = 12, col = 2, add = TRUE, lwd = 1.5)
Smith <- fitextcoeff(data, coord, estim = "Smith", loess = FALSE, ylim = c(1, 2))
plot(ext.coeff, from = 0, to = 12, col = 2, add = TRUE, lwd = 1.5)


###################################################
### code chunk number 21: SpatialExtremesGuide.Rnw:1424-1433 (eval = FALSE)
###################################################
## n.site <- 40
## n.obs <- 100
## coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
## data <- rmaxstab(n.obs, coord, "whitmat", nugget = 0, range = 1,
##                  smooth = 1)
## 
## par(mfrow=c(1,2))
## fitextcoeff(data, coord, loess = FALSE)
## fitextcoeff(data, coord, estim = "Smith", loess = FALSE)


###################################################
### code chunk number 22: madogram
###################################################
cov.fun1 <- covariance(nugget = 0, sill = 1, range = 1, smooth = 1,
                       cov.mod = "whitmat", plot = FALSE)
ext.coeff1 <- function(h)
  1 + sqrt((1 - cov.fun1(h))/2)
mado1 <- function(h)
  log(ext.coeff1(h))

cov.fun2 <- covariance(nugget = 0, sill = 1, range = 1.5, smooth = 1,
                       cov.mod = "powexp", plot = FALSE)
ext.coeff2 <- function(h)
  1 + sqrt((1 - cov.fun2(h))/2)
mado2 <- function(h)
  log(ext.coeff2(h))

n.site <- 40
n.obs <- 100
coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
data <- rmaxstab(n.obs, coord, "whitmat", nugget = 0, range = 1, smooth = 1)

par(mfrow=c(1,2))
madogram(data, coord, which = "mado", ylim = c(0, log(2)))
plot(mado1, from = 0, to = 12, add = TRUE, col = 2, lwd = 1.5)
madogram(data, coord, which = "mado", ylim = c(0, log(2)),
         n.bins = 100)
plot(mado1, from = 0, to = 12, add = TRUE, col = 2, lwd = 1.5)


###################################################
### code chunk number 23: SpatialExtremesGuide.Rnw:1519-1520
###################################################
getOption("SweaveHooks")[["fig"]]()
cov.fun1 <- covariance(nugget = 0, sill = 1, range = 1, smooth = 1,
                       cov.mod = "whitmat", plot = FALSE)
ext.coeff1 <- function(h)
  1 + sqrt((1 - cov.fun1(h))/2)
mado1 <- function(h)
  log(ext.coeff1(h))

cov.fun2 <- covariance(nugget = 0, sill = 1, range = 1.5, smooth = 1,
                       cov.mod = "powexp", plot = FALSE)
ext.coeff2 <- function(h)
  1 + sqrt((1 - cov.fun2(h))/2)
mado2 <- function(h)
  log(ext.coeff2(h))

n.site <- 40
n.obs <- 100
coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
data <- rmaxstab(n.obs, coord, "whitmat", nugget = 0, range = 1, smooth = 1)

par(mfrow=c(1,2))
madogram(data, coord, which = "mado", ylim = c(0, log(2)))
plot(mado1, from = 0, to = 12, add = TRUE, col = 2, lwd = 1.5)
madogram(data, coord, which = "mado", ylim = c(0, log(2)),
         n.bins = 100)
plot(mado1, from = 0, to = 12, add = TRUE, col = 2, lwd = 1.5)


###################################################
### code chunk number 24: SpatialExtremesGuide.Rnw:1552-1560 (eval = FALSE)
###################################################
## n.site <- 50
## n.obs <- 100
## coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
## data <- rmaxstab(n.obs, coord, "whitmat", nugget = 0, range = 1, smooth = 1)
## 
## par(mfrow=c(1,2))
## madogram(data, coord, which = "mado")
## madogram(data, coord, which = "mado", n.bins = 100)


###################################################
### code chunk number 25: madogramExtCoeff
###################################################
n.site <- 40
n.obs <- 100
coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
data <- rmaxstab(n.obs, coord, "whitmat", nugget = 0, range = 1, smooth = 1)

cov.fun <- covariance(nugget = 0, sill = 1, range = 1, smooth = 1, cov.mod = "whitmat",
                      plot = FALSE)
ext.coeff <- function(h)
  1 + sqrt((1 - cov.fun(h))/2)
mado <- function(h)
  log(ext.coeff(h))

par(mfrow=c(1,2))
madogram(data, coord, which = "mado", ylim = c(0, log(2)))
plot(mado, from = 0, to = 12, col = 2, lwd = 1.5, add = TRUE)
madogram(data, coord, which = "ext")
plot(ext.coeff, from = 0, to = 12, col = 2, lwd = 1.5, add = TRUE)


###################################################
### code chunk number 26: SpatialExtremesGuide.Rnw:1597-1598
###################################################
getOption("SweaveHooks")[["fig"]]()
n.site <- 40
n.obs <- 100
coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
data <- rmaxstab(n.obs, coord, "whitmat", nugget = 0, range = 1, smooth = 1)

cov.fun <- covariance(nugget = 0, sill = 1, range = 1, smooth = 1, cov.mod = "whitmat",
                      plot = FALSE)
ext.coeff <- function(h)
  1 + sqrt((1 - cov.fun(h))/2)
mado <- function(h)
  log(ext.coeff(h))

par(mfrow=c(1,2))
madogram(data, coord, which = "mado", ylim = c(0, log(2)))
plot(mado, from = 0, to = 12, col = 2, lwd = 1.5, add = TRUE)
madogram(data, coord, which = "ext")
plot(ext.coeff, from = 0, to = 12, col = 2, lwd = 1.5, add = TRUE)


###################################################
### code chunk number 27: SpatialExtremesGuide.Rnw:1615-1621 (eval = FALSE)
###################################################
## n.site <- 40
## n.obs <- 100
## coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
## data <- rmaxstab(n.obs, coord, "whitmat", nugget = 0, range = 1, smooth = 1)
## 
## madogram(data, coord)


###################################################
### code chunk number 28: F-madogram
###################################################
n.site <- 40
n.obs <- 100
coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
data <- rmaxstab(n.obs, coord, "gauss", cov11 = 1, cov12 = 0, cov22 = 1)

ext.coeff <- function(h)
  2 * pnorm(h/2)
Fmado <- function(h)
  0.5 * (ext.coeff(h) - 1) / (ext.coeff(h) + 1)

par(mfrow=c(1,2))
fmadogram(data, coord, which = "mado", ylim = c(0, 1/5))
plot(Fmado, from = 0, to = 13, add = TRUE, col = 2, lwd = 1.5)
fmadogram(data, coord, which = "ext")
plot(ext.coeff, from = 0, to = 13, add = TRUE, col = 2, lwd = 1.5)


###################################################
### code chunk number 29: SpatialExtremesGuide.Rnw:1660-1661
###################################################
getOption("SweaveHooks")[["fig"]]()
n.site <- 40
n.obs <- 100
coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
data <- rmaxstab(n.obs, coord, "gauss", cov11 = 1, cov12 = 0, cov22 = 1)

ext.coeff <- function(h)
  2 * pnorm(h/2)
Fmado <- function(h)
  0.5 * (ext.coeff(h) - 1) / (ext.coeff(h) + 1)

par(mfrow=c(1,2))
fmadogram(data, coord, which = "mado", ylim = c(0, 1/5))
plot(Fmado, from = 0, to = 13, add = TRUE, col = 2, lwd = 1.5)
fmadogram(data, coord, which = "ext")
plot(ext.coeff, from = 0, to = 13, add = TRUE, col = 2, lwd = 1.5)


###################################################
### code chunk number 30: SpatialExtremesGuide.Rnw:1744-1751 (eval = FALSE)
###################################################
## n.site <- 40
## n.obs <- 100
## coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
## data <- rmaxstab(n.obs, coord, "gauss", cov11 = 1, cov12 = 0, cov22 = 1)
## 
## par(mfrow=c(1,2))
## fmadogram(data, coord)


###################################################
### code chunk number 31: lambda-madogram
###################################################
n.site <- 40
n.obs <- 100
coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
data1 <- rmaxstab(n.obs, coord, "powexp", nugget = 0, range = 1, smooth = 1)
data2 <- rmaxstab(n.obs, coord, "cauchy", nugget = 0, range = 1, smooth = 1)

par(mfrow=c(1,2), pty = "s")
lmadogram(data1, coord, n.bins = 60)
lmadogram(data2, coord, n.bins = 60)


###################################################
### code chunk number 32: SpatialExtremesGuide.Rnw:1899-1900
###################################################
getOption("SweaveHooks")[["fig"]]()
n.site <- 40
n.obs <- 100
coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
data1 <- rmaxstab(n.obs, coord, "powexp", nugget = 0, range = 1, smooth = 1)
data2 <- rmaxstab(n.obs, coord, "cauchy", nugget = 0, range = 1, smooth = 1)

par(mfrow=c(1,2), pty = "s")
lmadogram(data1, coord, n.bins = 60)
lmadogram(data2, coord, n.bins = 60)


###################################################
### code chunk number 33: SpatialExtremesGuide.Rnw:1914-1923 (eval = FALSE)
###################################################
## n.site <- 40
## n.obs <- 100
## coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
## data1 <- rmaxstab(n.obs, coord, "powexp", nugget = 0, range = 1, smooth = 1)
## data2 <- rmaxstab(n.obs, coord, "cauchy", nugget = 0, range = 1, smooth = 1)
## 
## par(mfrow=c(1,2), pty = "s")
## lmadogram(data1, coord, n.bins = 60)
## lmadogram(data2, coord, n.bins = 60)


###################################################
### code chunk number 34: SpatialExtremesGuide.Rnw:1977-1983
###################################################
n.site <- 40
n.obs <- 80
coord <- matrix(runif(2*n.site, 0, 10), ncol = 2)
data <- rmaxstab(n.obs, coord, "gauss", cov11 = 9/8, cov12 = 1, cov22 = 9/8)

fitcovmat(data, coord, marge = "emp")


###################################################
### code chunk number 35: SpatialExtremesGuide.Rnw:2118-2125
###################################################
n.obs <- 80
n.site <- 40
set.seed(12)
coord <- matrix(runif(2*n.site, 0, 10), ncol = 2)
data <- rmaxstab(n.obs, coord, "whitmat", nugget = 0.8, range = 3, smooth = 1.2)

fitmaxstab(data, coord, cov.mod = "whitmat")


###################################################
### code chunk number 36: SpatialExtremesGuide.Rnw:2146-2149 (eval = FALSE)
###################################################
## fitmaxstab(data, coord, cov.mod = "whitmat", smooth = 1.2)
## fitmaxstab(data, coord, cov.mod = "whitmat", nugget = 0)
## fitmaxstab(data, coord, cov.mod = "whitmat", range = 3)


###################################################
### code chunk number 37: SpatialExtremesGuide.Rnw:2155-2157 (eval = FALSE)
###################################################
## fitmaxstab(data, coord, cov.mod = "cauchy")
## fitmaxstab(data, coord, cov.mod = "powexp")


###################################################
### code chunk number 38: SpatialExtremesGuide.Rnw:2161-2162 (eval = FALSE)
###################################################
## fitmaxstab(data, coord, cov.mod = "gauss")


###################################################
### code chunk number 39: SpatialExtremesGuide.Rnw:2168-2169 (eval = FALSE)
###################################################
## fitmaxstab(data, coord, cov.mod = "gauss", cov12 = 0, method = "BFGS")


###################################################
### code chunk number 40: SpatialExtremesGuide.Rnw:2341-2342
###################################################
set.seed(1)


###################################################
### code chunk number 41: SpatialExtremesGuide.Rnw:2345-2353
###################################################
n.obs <- 80
n.site <- 40
coord <- matrix(runif(2*n.site, 0, 10), ncol = 2)
data <- rmaxstab(n.obs, coord, "cauchy", nugget = 0.2, range = 3, smooth = 1.2)

M0 <- fitmaxstab(data, coord, cov.mod = "powexp")
M1 <- fitmaxstab(data, coord, cov.mod = "cauchy")
TIC(M0, M1)


###################################################
### code chunk number 42: SpatialExtremesGuide.Rnw:2419-2420
###################################################
set.seed(7)


###################################################
### code chunk number 43: SpatialExtremesGuide.Rnw:2423-2431
###################################################
n.obs <- 50
n.site <- 30
coord <- matrix(rnorm(2*n.site, sd = sqrt(.2)), ncol = 2)
data <- rmaxstab(n.obs, coord, "gauss", cov11 = 100, cov12 = 25, cov22 = 220)

M0 <- fitmaxstab(data, coord, "gauss", cov11 = 100)
M1 <- fitmaxstab(data, coord, "gauss")
anova(M0, M1)


###################################################
### code chunk number 44: SpatialExtremesGuide.Rnw:2509-2510
###################################################
anova(M0, M1, method = "CB")


###################################################
### code chunk number 45: SpatialExtremesGuide.Rnw:2545-2548
###################################################
x <- c(2.2975896, 1.6448808, 1.3323833, -0.4464904, 2.2737603, -0.2581876,
       9.5184398, -0.5899699, 0.4974283, -0.8152157)
z <- gev2frech(x, 1, 2, .2)


###################################################
### code chunk number 46: SpatialExtremesGuide.Rnw:2559-2560
###################################################
frech2gev(z, 1, 2, .2)


###################################################
### code chunk number 47: SpatialExtremesGuide.Rnw:2617-2618 (eval = FALSE)
###################################################
## fitmaxstab(data, coord, "gauss", fit.marge = TRUE)


###################################################
### code chunk number 48: knotsEffect
###################################################
set.seed(12)
x <- runif(100)
fun <- function(x) sin(3 * pi * x)
y <- fun(x) + rnorm(100, 0, 0.15)
knots1 <- quantile(x, prob = 1:2 / 3)
knots2 <- quantile(x, prob = 1:10 / 11)
knots3 <- quantile(x, prob = 1:50 / 51)
M0 <- rbpspline(y, x, knots = knots1, degree = 3, penalty = 0)
M1 <- rbpspline(y, x, knots = knots2, degree = 3, penalty = 0)
M2 <- rbpspline(y, x, knots = knots3, degree = 3, penalty = 0)

par(mfrow=c(1,3))
plot(x, y, col = "lightgrey")
rug(knots1)
lines(M0)
plot(x, y, col = "lightgrey")
rug(knots2)
lines(M1, col = 2)
plot(x, y, col = "lightgrey")
rug(knots3)
lines(M2, col = 3)


###################################################
### code chunk number 49: SpatialExtremesGuide.Rnw:2852-2853
###################################################
getOption("SweaveHooks")[["fig"]]()
set.seed(12)
x <- runif(100)
fun <- function(x) sin(3 * pi * x)
y <- fun(x) + rnorm(100, 0, 0.15)
knots1 <- quantile(x, prob = 1:2 / 3)
knots2 <- quantile(x, prob = 1:10 / 11)
knots3 <- quantile(x, prob = 1:50 / 51)
M0 <- rbpspline(y, x, knots = knots1, degree = 3, penalty = 0)
M1 <- rbpspline(y, x, knots = knots2, degree = 3, penalty = 0)
M2 <- rbpspline(y, x, knots = knots3, degree = 3, penalty = 0)

par(mfrow=c(1,3))
plot(x, y, col = "lightgrey")
rug(knots1)
lines(M0)
plot(x, y, col = "lightgrey")
rug(knots2)
lines(M1, col = 2)
plot(x, y, col = "lightgrey")
rug(knots3)
lines(M2, col = 3)


###################################################
### code chunk number 50: SpatialExtremesGuide.Rnw:2869-2870 (eval = FALSE)
###################################################
## set.seed(12)
## x <- runif(100)
## fun <- function(x) sin(3 * pi * x)
## y <- fun(x) + rnorm(100, 0, 0.15)
## knots1 <- quantile(x, prob = 1:2 / 3)
## knots2 <- quantile(x, prob = 1:10 / 11)
## knots3 <- quantile(x, prob = 1:50 / 51)
## M0 <- rbpspline(y, x, knots = knots1, degree = 3, penalty = 0)
## M1 <- rbpspline(y, x, knots = knots2, degree = 3, penalty = 0)
## M2 <- rbpspline(y, x, knots = knots3, degree = 3, penalty = 0)
## 
## par(mfrow=c(1,3))
## plot(x, y, col = "lightgrey")
## rug(knots1)
## lines(M0)
## plot(x, y, col = "lightgrey")
## rug(knots2)
## lines(M1, col = 2)
## plot(x, y, col = "lightgrey")
## rug(knots3)
## lines(M2, col = 3)


###################################################
### code chunk number 51: smoothingParamEffect
###################################################
M0 <- rbpspline(y, x, knots = knots3, degree = 3, penalty = 0)
M1 <- rbpspline(y, x, knots = knots3, degree = 3, penalty = 0.1)
M2 <- rbpspline(y, x, knots = knots3, degree = 3, penalty = 10)

par(mfrow=c(1,3))
plot(x, y, col = "lightgrey")
lines(M0)
plot(x, y, col = "lightgrey")
lines(M1, col = 2)
plot(x, y, col = "lightgrey")
lines(M2, col = 3)


###################################################
### code chunk number 52: SpatialExtremesGuide.Rnw:2895-2896
###################################################
getOption("SweaveHooks")[["fig"]]()
M0 <- rbpspline(y, x, knots = knots3, degree = 3, penalty = 0)
M1 <- rbpspline(y, x, knots = knots3, degree = 3, penalty = 0.1)
M2 <- rbpspline(y, x, knots = knots3, degree = 3, penalty = 10)

par(mfrow=c(1,3))
plot(x, y, col = "lightgrey")
lines(M0)
plot(x, y, col = "lightgrey")
lines(M1, col = 2)
plot(x, y, col = "lightgrey")
lines(M2, col = 3)


###################################################
### code chunk number 53: SpatialExtremesGuide.Rnw:2928-2929 (eval = FALSE)
###################################################
## M0 <- rbpspline(y, x, knots = knots3, degree = 3, penalty = 0)
## M1 <- rbpspline(y, x, knots = knots3, degree = 3, penalty = 0.1)
## M2 <- rbpspline(y, x, knots = knots3, degree = 3, penalty = 10)
## 
## par(mfrow=c(1,3))
## plot(x, y, col = "lightgrey")
## lines(M0)
## plot(x, y, col = "lightgrey")
## lines(M1, col = 2)
## plot(x, y, col = "lightgrey")
## lines(M2, col = 3)


###################################################
### code chunk number 54: CVandGCV
###################################################
par(mfrow=c(1,3))
lambda.cv <- cv(y, x, knots = knots3, degree = 3)$penalty
abline(v = lambda.cv, lty = 2)
lambda.gcv <- gcv(y, x, knots = knots3, degree = 3)$penalty
abline(v = lambda.gcv, lty = 2)
cv.fit <- rbpspline(y, x, knots3, degree = 3, penalty = "cv")
gcv.fit <- rbpspline(y, x, knots3, degree = 3, penalty = "gcv")
plot(x, y, col = "lightgrey")
lines(cv.fit, col = 2)
lines(gcv.fit, col = 3)


###################################################
### code chunk number 55: SpatialExtremesGuide.Rnw:3020-3021
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(1,3))
lambda.cv <- cv(y, x, knots = knots3, degree = 3)$penalty
abline(v = lambda.cv, lty = 2)
lambda.gcv <- gcv(y, x, knots = knots3, degree = 3)$penalty
abline(v = lambda.gcv, lty = 2)
cv.fit <- rbpspline(y, x, knots3, degree = 3, penalty = "cv")
gcv.fit <- rbpspline(y, x, knots3, degree = 3, penalty = "gcv")
plot(x, y, col = "lightgrey")
lines(cv.fit, col = 2)
lines(gcv.fit, col = 3)


###################################################
### code chunk number 56: SpatialExtremesGuide.Rnw:3036-3037 (eval = FALSE)
###################################################
## par(mfrow=c(1,3))
## lambda.cv <- cv(y, x, knots = knots3, degree = 3)$penalty
## abline(v = lambda.cv, lty = 2)
## lambda.gcv <- gcv(y, x, knots = knots3, degree = 3)$penalty
## abline(v = lambda.gcv, lty = 2)
## cv.fit <- rbpspline(y, x, knots3, degree = 3, penalty = "cv")
## gcv.fit <- rbpspline(y, x, knots3, degree = 3, penalty = "gcv")
## plot(x, y, col = "lightgrey")
## lines(cv.fit, col = 2)
## lines(gcv.fit, col = 3)


###################################################
### code chunk number 57: SpatialExtremesGuide.Rnw:3113-3114
###################################################
set.seed(15)


###################################################
### code chunk number 58: SpatialExtremesGuide.Rnw:3117-3130
###################################################
n.site <- 20
n.obs <- 50
coord <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(coord) <- c("lon", "lat")

data <- rmaxstab(n.obs, coord, "gauss", cov11 = 100, cov12 = 25, cov22 = 220)

param.loc <- -10 + 2 * coord[,2]
param.scale <- 5 + 2 * coord[,1] + coord[,2]^2
param.shape <- rep(0.2, n.site)

for (i in 1:n.site)
  data[,i] <- frech2gev(data[,i], param.loc[i], param.scale[i], param.shape[i])


###################################################
### code chunk number 59: SpatialExtremesGuide.Rnw:3135-3142
###################################################
loc.form <- y ~ lat
scale.form <- y ~ lon + I(lat^2)
shape.form <- y ~ 1
shape.form2 <- y ~ lon
M1 <- fitspatgev(data, coord, loc.form, scale.form, shape.form)
M2 <- fitspatgev(data, coord, loc.form, scale.form, shape.form2)
M1


###################################################
### code chunk number 60: SpatialExtremesGuide.Rnw:3149-3151
###################################################
anova(M1, M2)
TIC(M1, M2)


