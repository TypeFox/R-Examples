
library(gridGraphics)

plot.data.frame1 <- function() {
    plot(OrchardSprays[1], method = "jitter")
}

plot.data.frame2 <- function() {
    plot(OrchardSprays[c(4,1)])
}

plot.data.frame3 <- function() {
    plot(OrchardSprays)
}

plot.data.frame4 <- function() {
    plot(iris)
}

plot.data.frame5 <- function() {
    plot(iris[5:4])
}

plot.data.frame6 <- function() {
    plot(women)
}

plotdiff(expression(plot.data.frame1()), "plot.data.frame-1")
plotdiff(expression(plot.data.frame2()), "plot.data.frame-2", antialias=FALSE)
plotdiff(expression(plot.data.frame3()), "plot.data.frame-3",
         width=10, height=10)
plotdiff(expression(plot.data.frame4()), "plot.data.frame-4",
         width=12, height=12)
plotdiff(expression(plot.data.frame5()), "plot.data.frame-5")
plotdiff(expression(plot.data.frame6()), "plot.data.frame-6")

plotdiffResult()
