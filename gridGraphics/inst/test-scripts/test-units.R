
library(gridGraphics)

units1 <- function() {
    ## plot labels offset 0.12 inches to the right
    ## of plotted symbols in a plot
    with(mtcars, {
        plot(mpg, disp, pch = 19, main = "Motor Trend Cars")
        text(mpg + xinch(0.12), disp, row.names(mtcars),
             adj = 0, cex = .7, col = "blue")
    })
}

plotdiff(expression(units1()), "units-1")

plotdiffResult()
