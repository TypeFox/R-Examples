## =============================================================================
## Errors for the Euler method
## Figure 2.1 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================
require(shape)

# ------------------------------------------------------------------------------
# Figure 1. The euler method and the LTE and GTE
# ------------------------------------------------------------------------------
# main settings
windows(width = 8, height = 5)
par(mar = c(4.1, 5.1, 4.1, 2.1), las = 1) # las to have horizontal labels
par(mfrow = c(1, 2))

# function and derivates
fun <- function (x) 2*x^2
der <- function (x) 4*x    # 1st derivate
der2<- function (x) 4      # 2nd derivate
der3<- function (x) 12     # 3rd derivate

# points
x1  <- 1
hh  <- 1

x2  <- x1 + hh
y1  <- fun(x1)
y2  <- fun(x2)
cex.axis <- 1.0

# main plot function
plotf <- function()  {
 curve(fun(x), 1, 2.5, axes = FALSE, xlab = "", ylab = "", lwd = 2,
       ylim = c(0, 10), xlim = c(0.75, 2.5))

 axis(1, labels = FALSE, at = c(0, 3))
 axis(1, at = c(1, 2), 
      labels = c(expression(x[0]), expression(x[1] == x[0]+h)),
      font = 2, cex.axis = cex.axis)

 axis(2, labels=FALSE, at = c(-5, 35))
 axis(2, at = c(y1, ye, y2),
      labels = c(expression(y(x[0])==y[0]), "", expression(y(x[0]+h))),
      font = 2, cex.axis = cex.axis)

 points(x1, y1, pch = 16, col = "black", cex = 2)
 points(x2, y2, pch = 16, col = "black", cex = 2)
}

# 1-st order expansion
ye <- fun(x1) + hh*der(x1)

plotf()
mtext(side = 2, at = ye,expression(y[1]), line = 1, cex=cex.axis)

hhs <- seq(from = 0, to = hh, length.out = 100) 
xx  <- x1+hhs
lines(xx, fun(x1) + hhs*der(x1))

points(x2, ye, pch = 16, col = "darkgrey", cex = 2)
segments(x1, y1, x2, y1, lty = 2, col = "darkgrey")
segments(x2, y1, x2, ye, lty = 2, col = "darkgrey")

text(2.1, 4, expression(h*y[0] ^"'"), cex = cex.axis)
text(1.55, 3.9, "y '", cex = cex.axis)
text(1.55, 1.9, "h", cex = cex.axis)
text(2.1, 7, "LTE", cex = cex.axis)
 
dy <- 0.15
Arrows(x2, ye + dy, x2, y2 - dy, code = 1,
       arr.adj = 1, arr.type = "triangle",
       col = "grey", arr.length = 0.2)
Arrows(x2, ye + dy, x2, y2 - dy, code = 2,
       arr.adj = 1, arr.type = "triangle",
       col = "grey", arr.length = 0.2)

writelabel("A")

# the GTE

# true curve and points
pm <- 4
curve(fun(x), 1, pm, axes = FALSE, xlab = "", 
      ylab = "", lwd = 2, xlim = c(0.75, pm*1.1))
axis(1, labels = FALSE, at = 0:pm)
xlabels <- c(expression(x[0]), expression(x[1]), expression(x[2]),
             expression(x[3]), expression(x[4]))

axis(1, at = 1:pm, labels = xlabels[1:pm], font = 2, cex.axis = cex.axis)
xi <- 1:pm
yi <- fun(xi)
axis(2, labels = FALSE, at = c(0, yi))

ylabels <- c(expression(y(x[0])), expression(y(x[1])), expression(y(x[2])),
             expression(y(x[3])), expression(y(x[4]))) 
axis(2, at = yi, labels = ylabels[1:pm], font = 2, cex.axis = cex.axis)
points(xi, yi, pch = 16, col = "black", cex = 1)

yy <- vector(length = pm) 
xx <- vector(length = pm) 
xx[1] <- x1
yy[1] <- fun(x1)
h <- 1
for (i in 2:pm) {
    xx[i] <- xx[i-1] + h
    yy[i] <- yy[i-1] + h*der(xx[i-1])
}   
lines(xx, yy)
points(xx, yy, pch = 16, col = "grey")

segments(xx[2], yy[2], xx[2], fun(xx[2]), 
         col = "darkgrey", lty = 2, lwd = 2)
dx <- 0.1
# text(xx[2]+dx,0.5*(yy[2]+fun(xx[2])),"LTE",cex=0.75)

segments(xx[pm], yy[pm], xx[pm], fun(xx[pm]), 
         col = "darkgrey", lty = 2, lwd = 2)
dx <- 0.1
text(xx[pm] + dx, 0.5*(yy[pm]+fun(xx[pm])), "GTE", adj = 0)

writelabel("B")

