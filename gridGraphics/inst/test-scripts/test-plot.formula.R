
library(gridGraphics)

plot.formula1 <- function() {
    par(mfrow = c(2,1))
    plot(Ozone ~ Wind, data = airquality, pch = as.character(Month))
    plot(Ozone ~ Wind, data = airquality, pch = as.character(Month),
         subset = Month != 7)
}

plot.formula2 <- function() {
    ## text.formula() can be very natural:
    wb <- within(warpbreaks, {
        time <- seq_along(breaks); W.T <- wool:tension })
    plot(breaks ~ time, data = wb, type = "b")
}

plotdiff(expression(plot.formula1()), "plot.formula-1")
plotdiff(expression(plot.formula2()), "plot.formula-2")

plotdiffResult()
