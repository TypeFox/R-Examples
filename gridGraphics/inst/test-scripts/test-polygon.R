
library(gridGraphics)

polygon1 <- function() {
    x <- c(1:9, 8:1)
    y <- c(1, 2*(5:3), 2, -1, 17, 9, 8, 2:9)
    par(mfcol = c(3, 1))
    for(xpd in c(FALSE, TRUE, NA)) {
        plot(1:10, main = paste("xpd =", xpd))
        box("figure", col = "pink", lwd = 3)
        polygon(x, y, xpd = xpd, col = "orange", lty = 2, lwd = 2,
                border = "red")
    }
}

polygon2 <- function() {
    n <- 100
    xx <- c(0:n, n:0)
    set.seed(1)
    yy <- c(c(0, cumsum(stats::rnorm(n))), rev(c(0, cumsum(stats::rnorm(n)))))
    plot   (xx, yy, type = "n", xlab = "Time", ylab = "Distance")
    polygon(xx, yy, col = "gray", border = "red")
    title("Distance Between Brownian Motions")
}

polygon3 <- function() {
    # Multiple polygons from NA values
    # and recycling of col, border, and lty
    par(mfrow = c(2, 1), mar=c(5, 4, 2, 2))
    plot(c(1, 9), 1:2, type = "n")
    polygon(1:9, c(2,1,2,1,1,2,1,2,1),
            col = c("red", "blue"),
            border = c("green", "yellow"),
            lwd = 3, lty = c("dashed", "solid"))
    plot(c(1, 9), 1:2, type = "n")
    polygon(1:9, c(2,1,2,1,NA,2,1,2,1),
            col = c("red", "blue"),
            border = c("green", "yellow"),
            lwd = 3, lty = c("dashed", "solid"))
}

polygon4 <- function() {
    # Line-shaded polygons
    plot(c(1, 9), 1:2, type = "n")
    polygon(1:9, c(2,1,2,1,NA,2,1,2,1),
            density = c(10, 20), angle = c(-45, 45))
}

plotdiff(expression(polygon1()), "polygon-1")
plotdiff(expression(polygon2()), "polygon-2")
plotdiff(expression(polygon3()), "polygon-3")
plotdiff(expression(polygon4()), "polygon-4")

plotdiffResult()
