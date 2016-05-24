
library(gridGraphics)

require(stats)

plot.design1 <- function() {
    plot.design(warpbreaks)  # automatic for data frame with one numeric var.
}

plot.design2 <- function() {
    Form <- breaks ~ wool + tension
    summary(fm1 <- aov(Form, data = warpbreaks))
    plot.design(       Form, data = warpbreaks, col = 2)  # same as above
}

plot.design3 <- function() {
    ## More than one y :
    plot.design(esoph) ## two plots; if interactive you are "ask"ed
}

plot.design4 <- function() {
    ## or rather, compare mean and median:
    par(mfcol = 1:2)
    plot.design(ncases/ncontrols ~ ., data = esoph, ylim = c(0, 0.8))
    plot.design(ncases/ncontrols ~ ., data = esoph, ylim = c(0, 0.8),
                fun = median)
}

plotdiff(expression(plot.design1()), "plot.design-1")
plotdiff(expression(plot.design2()), "plot.design-2")
plotdiff(expression(plot.design3()), "plot.design-3")
plotdiff(expression(plot.design4()), "plot.design-4", width=10)

plotdiffResult()
