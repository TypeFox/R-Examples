
library(gridGraphics)

require(stats)

plot1 <- function() {
    plot(cars)
    lines(lowess(cars))
}

plot2 <- function() {
    plot(sin, -pi, 2*pi) # see ?plot.function
}

plot3 <- function() {
    ## Discrete Distribution Plot:
    plot(table(rpois(100, 5)), type = "h", col = "red", lwd = 10,
         main = "rpois(100, lambda = 5)")
}

plot4 <- function() {
    ## Simple quantiles/ECDF, see ecdf() {library(stats)} for a better one:
    plot(x <- sort(rnorm(47)), type = "s", main = "plot(x, type = \"s\")")
    points(x, cex = .5, col = "dark red")
}

plotdiff(expression(plot1()), "plot-1")
plotdiff(expression(plot2()), "plot-2")
plotdiff(expression(plot3()), "plot-3")
plotdiff(expression(plot4()), "plot-4")

plotdiffResult()
