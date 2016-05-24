### R code from vignette source 'Ch-Viz.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library("MVA")
set.seed(280875)
library("lattice")
lattice.options(default.theme =
    function()
        standard.theme("pdf", color = FALSE))

if (file.exists("deparse.R")) {
    if (!file.exists("figs")) dir.create("figs")
    source("deparse.R")
    options(prompt = "R> ", continue = "+  ", width = 64,
        digits = 4, show.signif.stars = FALSE, useFancyQuotes = FALSE)

    options(SweaveHooks = list(onefig =   function() {par(mfrow = c(1,1))},
                           twofig =   function() {par(mfrow = c(1,2))},                           
                           figtwo =   function() {par(mfrow = c(2,1))},                           
                           threefig = function() {par(mfrow = c(1,3))},
                           figthree = function() {par(mfrow = c(3,1))},
                           fourfig =  function() {par(mfrow = c(2,2))},
                           sixfig =   function() {par(mfrow = c(3,2))},
                           nomar = function() par("mai" = c(0, 0, 0, 0))))
}


###################################################
### code chunk number 2: ch:Viz:data
###################################################
measure <- 
structure(list(V1 = 1:20, V2 = c(34L, 37L, 38L, 36L, 38L, 43L,
40L, 38L, 40L, 41L, 36L, 36L, 34L, 33L, 36L, 37L, 34L, 36L, 38L,
35L), V3 = c(30L, 32L, 30L, 33L, 29L, 32L, 33L, 30L, 30L, 32L,
24L, 25L, 24L, 22L, 26L, 26L, 25L, 26L, 28L, 23L), V4 = c(32L,
37L, 36L, 39L, 33L, 38L, 42L, 40L, 37L, 39L, 35L, 37L, 37L, 34L,
38L, 37L, 38L, 37L, 40L, 35L)), .Names = c("V1", "V2", "V3",
"V4"), class = "data.frame", row.names = c(NA, -20L))
measure <- measure[,-1]
names(measure) <- c("chest", "waist", "hips")
measure$gender <- gl(2, 10)
levels(measure$gender) <- c("male", "female")

data("USairpollution", package = "HSAUR2")


###################################################
### code chunk number 3: ch:Viz:USairpollution:plot1mlab
###################################################
mlab <- "Manufacturing enterprises with 20 or more workers"
plab <- "Population size (1970 census) in thousands"


###################################################
### code chunk number 4: ch:Viz:USairpollution:plot1setup (eval = FALSE)
###################################################
## plot(popul ~ manu, data = USairpollution, 
##      xlab = mlab, ylab = plab)


###################################################
### code chunk number 5: ch:Viz:USairpollution:plot1
###################################################
plot(popul ~ manu, data = USairpollution, 
     xlab = mlab, ylab = plab)


###################################################
### code chunk number 6: ch:Viz:USairpollution:plot3setup
###################################################
layout(matrix(c(2, 0, 1, 3), nrow = 2, byrow = TRUE),
       widths = c(2, 1), heights = c(1, 2), respect = TRUE)
xlim <- with(USairpollution, range(manu)) * 1.1
plot(popul ~ manu, data = USairpollution, cex.lab = 0.9,
     xlab = mlab, ylab = plab, type = "n", xlim = xlim)
with(USairpollution, text(manu, popul, cex = 0.6,
     labels = abbreviate(row.names(USairpollution))))
with(USairpollution, hist(manu, main = "", xlim = xlim))
with(USairpollution, boxplot(popul))


###################################################
### code chunk number 7: ch:Viz:USairpollution:plot2
###################################################
plot(popul ~ manu, data = USairpollution, 
     xlab = mlab, ylab = plab)
rug(USairpollution$manu, side = 1)
rug(USairpollution$popul, side = 2)


###################################################
### code chunk number 8: ch:Viz:USairpollution:plot3
###################################################
layout(matrix(c(2, 0, 1, 3), nrow = 2, byrow = TRUE),
       widths = c(2, 1), heights = c(1, 2), respect = TRUE)
xlim <- with(USairpollution, range(manu)) * 1.1
plot(popul ~ manu, data = USairpollution, cex.lab = 0.9,
     xlab = mlab, ylab = plab, type = "n", xlim = xlim)
with(USairpollution, text(manu, popul, cex = 0.6,
     labels = abbreviate(row.names(USairpollution))))
with(USairpollution, hist(manu, main = "", xlim = xlim))
with(USairpollution, boxplot(popul))


###################################################
### code chunk number 9: ch:Viz:USairpollution:plot4
###################################################
outcity <- match(lab <- c("Chicago", "Detroit", 
    "Cleveland", "Philadelphia"), rownames(USairpollution))
x <- USairpollution[, c("manu", "popul")]
bvbox(x, mtitle = "", xlab = mlab, ylab = plab)
text(x$manu[outcity], x$popul[outcity], labels = lab,
     cex = 0.7, pos = c(2, 2, 4, 2, 2))


###################################################
### code chunk number 10: ch:Viz:USairpollution:cor
###################################################
with(USairpollution, cor(manu, popul))
outcity <- match(c("Chicago", "Detroit", 
                   "Cleveland", "Philadelphia"),
                 rownames(USairpollution))
with(USairpollution, cor(manu[-outcity], popul[-outcity]))


###################################################
### code chunk number 11: ch:Viz:USairpollution:chull
###################################################
(hull <- with(USairpollution, chull(manu, popul)))


###################################################
### code chunk number 12: ch:Viz:USairpollution:chullplot
###################################################
with(USairpollution, 
     plot(manu, popul, pch = 1, xlab = mlab, ylab = plab))
with(USairpollution, 
     polygon(manu[hull], popul[hull], density = 15, angle = 30))


###################################################
### code chunk number 13: ch:Viz:USairpollution:chullcor
###################################################
with(USairpollution, cor(manu[-hull],popul[-hull]))


###################################################
### code chunk number 14: ch:Viz:USairpollution:chiplot:setup (eval = FALSE)
###################################################
## with(USairpollution, plot(manu, popul, 
##                           xlab = mlab, ylab = plab, 
##                           cex.lab = 0.9))
## with(USairpollution, chiplot(manu, popul))


###################################################
### code chunk number 15: ch:Viz:USairpollution:chiplot
###################################################
with(USairpollution, plot(manu, popul, 
                          xlab = mlab, ylab = plab, 
                          cex.lab = 0.9))
with(USairpollution, chiplot(manu, popul))


###################################################
### code chunk number 16: ch:Viz:USairpollution:plot5
###################################################
ylim <- with(USairpollution, range(wind)) * c(0.95, 1)
plot(wind ~ temp, data = USairpollution, 
     xlab = "Average annual temperature (Fahrenheit)",
     ylab = "Average annual wind speed (m.p.h.)", pch = 10,
     ylim = ylim)
with(USairpollution, symbols(temp, wind, circles = SO2, 
                             inches = 0.5, add = TRUE))


###################################################
### code chunk number 17: ch:Viz:USairpollution:plot6
###################################################
plot(wind ~ temp, data = USairpollution,
     xlab = "Average annual temperature (Fahrenheit)",
     ylab = "Average annual wind speed (m.p.h.)", pch = 10,
     ylim = ylim)
with(USairpollution,
    stars(USairpollution[,-c(2,5)], locations = cbind(temp, wind),
          labels = NULL, add = TRUE, cex = 0.5))


###################################################
### code chunk number 18: ch:Viz:USairpollution:plot7
###################################################
stars(USairpollution, cex = 0.55)


###################################################
### code chunk number 19: ch:Viz:USairpollution:plot8
###################################################
pairs(USairpollution, pch = ".", cex = 1.5)


###################################################
### code chunk number 20: ch:Viz:USairpollution:plot9
###################################################
pairs(USairpollution, 
      panel = function (x, y, ...) {
          points(x, y, ...)
          abline(lm(y ~ x), col = "grey")
      }, pch = ".", cex = 1.5)


###################################################
### code chunk number 21: ch:Viz:USairpollution:cor
###################################################
round(cor(USairpollution), 4) 


###################################################
### code chunk number 22: ch:Viz-kernel-figs
###################################################
rec <- function(x) (abs(x) < 1) * 0.5
tri <- function(x) (abs(x) < 1) * (1 - abs(x))
gauss <- function(x) 1/sqrt(2*pi) * exp(-(x^2)/2)
x <- seq(from = -3, to = 3, by = 0.001)
plot(x, rec(x), type = "l", ylim = c(0,1), lty = 1, 
     ylab = expression(K(x)))
lines(x, tri(x), lty = 2)
lines(x, gauss(x), lty = 3)
legend("topleft", legend = c("Rectangular", "Triangular", 
       "Gaussian"), lty = 1:3, title = "kernel functions", 
       bty = "n")


###################################################
### code chunk number 23: ch:Viz-x-bumps-data
###################################################
x <- c(0, 1, 1.1, 1.5, 1.9, 2.8, 2.9, 3.5)
n <- length(x)


###################################################
### code chunk number 24: ch:Viz-x-bumps-gaussian
###################################################
xgrid <- seq(from = min(x) - 1, to = max(x) + 1, by = 0.01) 


###################################################
### code chunk number 25: ch:Viz-x-bumps-bumps
###################################################
h <- 0.4
bumps <- sapply(x, function(a) gauss((xgrid - a)/h)/(n * h))


###################################################
### code chunk number 26: ch:Viz-x-bumps-setup (eval = FALSE)
###################################################
## plot(xgrid, rowSums(bumps), ylab = expression(hat(f)(x)),
##      type = "l", xlab = "x", lwd = 2)
## rug(x, lwd = 2)
## out <- apply(bumps, 2, function(b) lines(xgrid, b))


###################################################
### code chunk number 27: ch:Viz-x-bumps
###################################################
plot(xgrid, rowSums(bumps), ylab = expression(hat(f)(x)),
     type = "l", xlab = "x", lwd = 2)
rug(x, lwd = 2)
out <- apply(bumps, 2, function(b) lines(xgrid, b))
par(op)


###################################################
### code chunk number 28: ch:Viz-epakernel-fig
###################################################
epa <- function(x, y) 
    ((x^2 + y^2) < 1) * 2/pi * (1 - x^2 - y^2)
x <- seq(from = -1.1, to = 1.1, by = 0.05)
epavals <- sapply(x, function(a) epa(a, x))
persp(x = x, y = x, z = epavals, xlab = "x", ylab = "y", 
      zlab = expression(K(x, y)), theta = -35, axes = TRUE, 
      box = TRUE)


###################################################
### code chunk number 29: ch:Viz:CYGOB1:tab
###################################################
toLatex(HSAURtable(CYGOB1), pcol = 3,
    caption = "Energy output and surface temperature of star cluster CYG OB1.",
    label = "ch:Viz:CYGOB1:tab")


###################################################
### code chunk number 30: ch:Viz:CYGOB1:plot1
###################################################
library("KernSmooth")
CYGOB1d <- bkde2D(CYGOB1, bandwidth = sapply(CYGOB1, dpik))
plot(CYGOB1, xlab = "log surface temperature",
             ylab = "log light intensity")
contour(x = CYGOB1d$x1, y = CYGOB1d$x2, 
        z = CYGOB1d$fhat, add = TRUE)


###################################################
### code chunk number 31: ch:Viz:CYGOB1:plot2
###################################################
persp(x = CYGOB1d$x1, y = CYGOB1d$x2, z = CYGOB1d$fhat,
      xlab = "log surface temperature",
      ylab = "log light intensity",
      zlab = "density")


###################################################
### code chunk number 32: ch:Viz:measure:plot1:setup (eval = FALSE)
###################################################
## panel.hist <- function(x, ...)
## {
##     usr <- par("usr"); on.exit(par(usr))
##     par(usr = c(usr[1:2], 0, 1.5) )
##     h <- hist(x, plot = FALSE)
##     breaks <- h$breaks; nB <- length(breaks)
##     y <- h$counts; y <- y/max(y)
##     rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
## }
## pairs(measure[, c("chest", "waist", "hips")],
##       diag.panel = panel.hist,
##       panel = function (x,y) {
##           data <- data.frame(cbind(x,y))
##           par(new = TRUE)
##           den <- bkde2D(data, bandwidth = sapply(data, dpik))
##           contour(x = den$x1, y = den$x2, 
##                   z = den$fhat, axes = FALSE)
##       })


###################################################
### code chunk number 33: ch:Viz:measure:plot1
###################################################
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}
pairs(measure[, c("chest", "waist", "hips")],
      diag.panel = panel.hist,
      panel = function (x,y) {
          data <- data.frame(cbind(x,y))
          par(new = TRUE)
          den <- bkde2D(data, bandwidth = sapply(data, dpik))
          contour(x = den$x1, y = den$x2, 
                  z = den$fhat, axes = FALSE)
      })


###################################################
### code chunk number 34: ch:Viz:measure:plot2
###################################################
library("scatterplot3d")
with(measure, scatterplot3d(chest, waist, hips,
     pch = (1:2)[gender], type = "h", angle = 55))


###################################################
### code chunk number 35: ch:Viz:USairpollution:plot10
###################################################
with(USairpollution, 
    scatterplot3d(temp, wind, SO2, type = "h",
                  angle = 55))


###################################################
### code chunk number 36: ch:Viz:USairpollution:plot11
###################################################
plot(xyplot(SO2 ~ temp| cut(wind, 2), data = USairpollution))


###################################################
### code chunk number 37: ch:Viz:USairpollution:plot12
###################################################
pollution <- with(USairpollution, equal.count(SO2,4))
plot(cloud(precip ~ temp * wind | pollution, panel.aspect = 0.9,
     data = USairpollution))


###################################################
### code chunk number 38: ch:Viz:quakes:plot1
###################################################
plot(xyplot(lat ~ long| cut(depth, 3), data = quakes, 
            layout = c(3, 1), xlab = "Longitude", 
            ylab = "Latitude"))


###################################################
### code chunk number 39: ch:Viz:quakes:plot2
###################################################
Magnitude <- with(quakes, equal.count(mag, 4))
depth.ord <- with(quakes, rev(order(depth)))
quakes.ordered <- quakes[depth.ord,]
depth.breaks <- with(quakes.ordered, do.breaks(range(depth),50))
quakes.ordered$color<-level.colors(quakes.ordered$depth,at=depth.breaks,
   col.regions=grey.colors)
plot(xyplot(lat ~ long | Magnitude, data = quakes.ordered,
       aspect = "iso", groups = color, cex = 2, col = "black",
       panel = function(x, y, groups, ..., subscripts) {
           fill <- groups[subscripts]
           panel.grid(h = -1, v = -1)
           panel.xyplot(x, y, pch = 21, fill = fill, ...)
       },
       legend =
       list(right =
            list(fun = draw.colorkey,
                 args = list(key = list(col = gray.colors,
                                        at = depth.breaks),
                             draw = FALSE))),
       xlab = "Longitude", ylab = "Latitude"))


###################################################
### code chunk number 40: ch:Viz:quakes:plot3
###################################################
plot(cloud(depth ~ lat * long | Magnitude, data = quakes,
      zlim = rev(range(quakes$depth)),
      screen = list(z = 105, x = -70), panel.aspect = 0.9,
      xlab = "Longitude", ylab = "Latitude", zlab = "Depth"))


###################################################
### code chunk number 41: ch:Viz:USairpollution:stalac
###################################################
stalac(USairpollution)


