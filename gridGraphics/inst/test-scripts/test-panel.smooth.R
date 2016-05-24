
library(gridGraphics)

panel.smooth1 <- function() {
    pairs(swiss, panel = panel.smooth, pch = ".")  # emphasize the smooths
}

panel.smooth2 <- function() {
    pairs(swiss, panel = panel.smooth, lwd = 2, cex = 1.5, col = "blue")  # hmm...
}

plotdiff(expression(panel.smooth1()), "panel.smooth-1", width=12, height=12)
plotdiff(expression(panel.smooth2()), "panel.smooth-2", width=12, height=12)

plotdiffResult()
