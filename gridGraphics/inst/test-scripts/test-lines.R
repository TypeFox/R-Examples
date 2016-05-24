
library(gridGraphics)

lines1 <- function() {
    plot(cars, main = "Stopping Distance versus Speed")
    lines(stats::lowess(cars))
}

plotdiff(expression(lines1()), "lines-1")

plotdiffResult()
