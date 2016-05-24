
library(gridGraphics)

require(stats)
require(grDevices)

sunflowerplot1 <- function() {
    ## 'number' is computed automatically:
    sunflowerplot(iris[, 3:4])
}

sunflowerplot2 <- function() {
    ## Imitating Chambers et al, p.109, closely:
    sunflowerplot(iris[, 3:4], cex = .2, cex.fact = 1, size = .035, seg.lwd = .8)
}

sunflowerplot3 <- function() {
    ## or
    sunflowerplot(Petal.Width ~ Petal.Length, data = iris,
                  cex = .2, cex.fact = 1, size = .035, seg.lwd = .8)
}

sunflowerplot4 <- function() {
    set.seed(1)
    sunflowerplot(x = sort(2*round(rnorm(100))), y = round(rnorm(100), 0),
                  main = "Sunflower Plot of Rounded N(0,1)")
}

## Similarly using a "xyTable" argument:
xyT <- xyTable(x = sort(2*round(rnorm(100))), y = round(rnorm(100), 0),
               digits = 3)
utils::str(xyT, vec.len = 20)

sunflowerplot5 <- function() {
    sunflowerplot(xyT, main = "2nd Sunflower Plot of Rounded N(0,1)")
}

sunflowerplot6 <- function() {
    ## A 'marked point process' {explicit 'number' argument}:
    set.seed(1)
    sunflowerplot(rnorm(100), rnorm(100), number = rpois(n = 100, lambda = 2),
                  main = "Sunflower plot (marked point process)",
                  rotate = TRUE, col = "blue4")
}

plotdiff(expression(sunflowerplot1()), "sunflowerplot-1")
plotdiff(expression(sunflowerplot2()), "sunflowerplot-2")
plotdiff(expression(sunflowerplot3()), "sunflowerplot-3")
plotdiff(expression(sunflowerplot4()), "sunflowerplot-4")
plotdiff(expression(sunflowerplot5()), "sunflowerplot-5")
plotdiff(expression(sunflowerplot6()), "sunflowerplot-6")

plotdiffResult()
