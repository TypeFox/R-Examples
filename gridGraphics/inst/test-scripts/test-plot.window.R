
require(stats)  # normally loaded

library(gridGraphics)

plot.window1 <- function() {
    loc <- cmdscale(eurodist)
    rx <- range(x <- loc[,1])
    ry <- range(y <- -loc[,2])
    plot(x, y, type = "n", asp = 1, xlab = "", ylab = "")
    abline(h = pretty(rx, 10), v = pretty(ry, 10), col = "lightgray")
    text(x, y, labels(eurodist), cex = 0.8)
}

plotdiff(expression(plot.window1()), "plot.window-1")

plotdiffResult()
