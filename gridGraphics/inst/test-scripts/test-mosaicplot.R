
library(gridGraphics)

require(stats)

mosaicplot1 <- function() {
    mosaicplot(Titanic, main = "Survival on the Titanic", color = TRUE)
}

mosaicplot2 <- function() {
    ## Formula interface for tabulated data:
    mosaicplot(~ Sex + Age + Survived, data = Titanic, color = TRUE)
}

mosaicplot3 <- function() {
    mosaicplot(HairEyeColor, shade = TRUE)
}

mosaicplot4 <- function() {
    mosaicplot(HairEyeColor, shade = TRUE, margin = list(1:2, 3))
}

mosaicplot5 <- function() {
    ## Formula interface for raw data: visualize cross-tabulation of numbers
    ## of gears and carburettors in Motor Trend car data.
    mosaicplot(~ gear + carb, data = mtcars, color = TRUE, las = 1)
}

mosaicplot6 <- function() {
    # color recycling
    mosaicplot(~ gear + carb, data = mtcars, color = 2:3, las = 1)
}

plotdiff(expression(mosaicplot1()), "mosaicplot-1")
plotdiff(expression(mosaicplot2()), "mosaicplot-2")
plotdiff(expression(mosaicplot3()), "mosaicplot-3")
plotdiff(expression(mosaicplot4()), "mosaicplot-4")
plotdiff(expression(mosaicplot5()), "mosaicplot-5")
plotdiff(expression(mosaicplot6()), "mosaicplot-6")

plotdiffResult()
