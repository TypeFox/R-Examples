
library(gridGraphics)

abline1 <- function() {
    ## Setup up coordinate system (with x == y aspect ratio):
    plot(c(-2,3), c(-1,5), type = "n", xlab = "x", ylab = "y", asp = 1)
    ## the x- and y-axis, and an integer grid
    abline(h = 0, v = 0, col = "gray60")
    text(1,0, "abline( h = 0 )", col = "gray60", adj = c(0, -.1))
    abline(h = -1:5, v = -2:3, col = "lightgray", lty = 3)
    abline(a = 1, b = 2, col = 2)
    text(1,3, "abline( 1, 2 )", col = 2, adj = c(-.1, -.1))
}

abline2 <- function() {
    ## Simple Regression Lines:
    require(stats)
    sale5 <- c(6, 4, 9, 7, 6, 12, 8, 10, 9, 13)
    plot(sale5)
    abline(lsfit(1:10, sale5))
    abline(lsfit(1:10, sale5, intercept = FALSE), col = 4) # less fitting
}

abline3 <- function() {
    z <- lm(dist ~ speed, data = cars)
    plot(cars)
    abline(z) # equivalent to abline(reg = z) or
    abline(coef = coef(z))

    ## trivial intercept model
    abline(mC <- lm(dist ~ 1, data = cars)) ## the same as
    abline(a = coef(mC), b = 0, col = "blue")
}

# Test 'untf' and log scales
abline4 <- function() {
    par(mfrow=c(2, 2), mar=c(5, 4, 2, 2))
    plot(1:10)
    abline(1, 1)
    plot(1:10, log="x")
    abline(1, 1, untf=TRUE)
    plot(1:10, log="y")
    abline(1, 1, untf=TRUE)
    plot(1:10, log="xy")
    abline(1, 1, untf=TRUE)
}

plotdiff(expression(abline1()), "abline-1")
plotdiff(expression(abline2()), "abline-2")
plotdiff(expression(abline3()), "abline-3")
plotdiff(expression(abline4()), "abline-4")

plotdiffResult()
