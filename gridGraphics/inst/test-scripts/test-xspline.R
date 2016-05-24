library(gridGraphics)

xsplineTest <- function(s, open = TRUE,
                        x = c(1,1,3,3)/4,
                        y = c(1,3,3,1)/4, ...) {
    plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "")
    points(x, y, pch = 19)
    xspline(x, y, s, open, ...)
    text(x+0.05*c(-1,-1,1,1), y+0.05*c(-1,1,1,-1), s)
}

xspline1 <- function() {
    par(mfrow = c(3,3), mar = rep(0,4), oma = c(0,0,2,0))
    xsplineTest(c(0, -1, -1, 0))
    xsplineTest(c(0, -1,  0, 0))
    xsplineTest(c(0, -1,  1, 0))
    xsplineTest(c(0,  0, -1, 0))
    xsplineTest(c(0,  0,  0, 0))
    xsplineTest(c(0,  0,  1, 0))
    xsplineTest(c(0,  1, -1, 0))
    xsplineTest(c(0,  1,  0, 0))
    xsplineTest(c(0,  1,  1, 0))
    title("Open X-splines", outer = TRUE)
}

xspline2 <- function() {
    par(mfrow = c(3,3), mar = rep(0,4), oma = c(0,0,2,0))
    xsplineTest(c(0, -1, -1, 0), FALSE, col = "grey80")
    xsplineTest(c(0, -1,  0, 0), FALSE, col = "grey80")
    xsplineTest(c(0, -1,  1, 0), FALSE, col = "grey80")
    xsplineTest(c(0,  0, -1, 0), FALSE, col = "grey80")
    xsplineTest(c(0,  0,  0, 0), FALSE, col = "grey80")
    xsplineTest(c(0,  0,  1, 0), FALSE, col = "grey80")
    xsplineTest(c(0,  1, -1, 0), FALSE, col = "grey80")
    xsplineTest(c(0,  1,  0, 0), FALSE, col = "grey80")
    xsplineTest(c(0,  1,  1, 0), FALSE, col = "grey80")
    title("Closed X-splines", outer = TRUE)
}

xspline3 <- function() {
    set.seed(1)
    x <- sort(stats::rnorm(5))
    y <- sort(stats::rnorm(5))
    plot(x, y, pch = 19)
    res <- xspline(x, y, 1, draw = FALSE)
    lines(res)
    ## the end points may be very close together,
    ## so use last few for direction
    nr <- length(res$x)
    arrows(res$x[1], res$y[1], res$x[4], res$y[4], code = 1, length = 0.1)
    arrows(res$x[nr-3], res$y[nr-3], res$x[nr], res$y[nr], code = 2,
           length = 0.1)
}

plotdiff(expression(xspline1()), "xspline-1")
plotdiff(expression(xspline2()), "xspline-2")
plotdiff(expression(xspline3()), "xspline-3")

plotdiffResult()

