
library(gridGraphics)

require(grDevices)

plot.factor1 <- function() {
    plot(weight ~ group, data = PlantGrowth)           # numeric vector ~ factor
}

plot.factor2 <- function() {
    plot(cut(weight, 2) ~ group, data = PlantGrowth)   # factor ~ factor
    ## passing "..." to spineplot() eventually:
}

plot.factor3 <- function() {
    plot(cut(weight, 3) ~ group, data = PlantGrowth,
         col = hcl(c(0, 120, 240), 50, 70))
}

plot.factor4 <- function() {
    plot(PlantGrowth$group, axes = FALSE, main = "no axes")  # extremely silly
}

plotdiff(expression(plot.factor1()), "plot.factor-1")
plotdiff(expression(plot.factor2()), "plot.factor-2")
plotdiff(expression(plot.factor3()), "plot.factor-3")
plotdiff(expression(plot.factor4()), "plot.factor-4")

plotdiffResult()
