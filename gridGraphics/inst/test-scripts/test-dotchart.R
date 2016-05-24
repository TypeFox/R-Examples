
library(gridGraphics)

dotchart1 <- function() {
    dotchart(VADeaths, main = "Death Rates in Virginia - 1940")
}

dotchart2 <- function() {
    par(xaxs = "i")  # 0 -- 100\%
    dotchart(t(VADeaths), xlim = c(0,100),
             main = "Death Rates in Virginia - 1940")
}

plotdiff(expression(dotchart1()), "dotchart-1")
plotdiff(expression(dotchart2()), "dotchart-2")

plotdiffResult()
