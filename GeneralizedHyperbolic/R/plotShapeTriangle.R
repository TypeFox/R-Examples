### Function to plot the shape triangle
plotShapeTriangle <- function(xgap = 0.025, ygap = 0.0625/2,
                              main = "Shape Triangle", ...){
  x <- 0
  y <- 0
  par(mar = c(5,5,4,2) + 0.1)
  xlim <- c(-1 - xgap, 1 + xgap)
  ylim <- c(-0.25 - ygap, 1.25 + ygap)
  plot(x, y, type = "n", xlim = xlim, ylim = ylim,
       xaxt = "n", yaxt = "n",
       xlab = expression(chi), ylab = "", las = 1,
       main = main, ...)
  axis(1, at = (0:11)*0.20 - 1.00, las = 1)
  axis(2, at = (0:6)*0.25 - 0.25, las = 1)
  mtext(expression(xi), side = 2, line = 4, las = 1)
  segments((0:11)*0.20 - 1, -0.25, (0:11)*0.20 - 1, 1.25, lty = 3)
  segments(-1, (0:6)*0.25 - 0.25, 1, (0:6)*0.25 - 0.25, lty = 3)
  segments(-1, 1, 1, 1, lty = 1, lwd = 1.3)
  segments(-1, 1, 0, 0, lty = 1, lwd = 1.3)
  segments(0, 0, 1, 1, lty = 1, lwd = 1.3)
}
