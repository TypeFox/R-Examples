
library(gridGraphics)

box1 <- function() {
    set.seed(1)
    plot(1:7, abs(stats::rnorm(7)), type = "h", axes = FALSE)
    axis(1, at = 1:7, labels = letters[1:7])
    box(lty = '1373', col = 'red')
}

plotdiff(expression(box1()), "box-1")

plotdiffResult()
