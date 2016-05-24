
library(gridGraphics)

library(KernSmooth)

## A largish data set
set.seed(1)
n <- 10000
x1  <- matrix(rnorm(n), ncol = 2)
x2  <- matrix(rnorm(n, mean = 3, sd = 1.5), ncol = 2)
x   <- rbind(x1, x2)

smoothScatter1 <- function() {
    par(mfrow = c(2, 2))
    smoothScatter(x, nrpoints = 0)
    smoothScatter(x)

    ## a different color scheme:
    Lab.palette <- colorRampPalette(c("blue", "orange", "red"), space = "Lab")
    smoothScatter(x, colramp = Lab.palette)

    ## somewhat similar, using identical smoothing computations,
    ## but considerably *less* efficient for really large data:
    plot(x, col = densCols(x), pch = 20)
}

smoothScatter2 <- function() {
## use with pairs:
    set.seed(1)
    y <- matrix(rnorm(40000), ncol = 4) + 3*rnorm(10000)
    y[, c(2,4)] <-  -y[, c(2,4)]
    pairs(y, panel = function(...) smoothScatter(..., nrpoints = 0, add = TRUE))
}

plotdiff(expression(smoothScatter1()), "smoothScatter-1", antialias=FALSE)
plotdiff(expression(smoothScatter2()), "smoothScatter-2",
         antialias=FALSE, width=10, height=10)

plotdiffResult()
