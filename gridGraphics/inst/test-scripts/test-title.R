
library(gridGraphics)

title1 <- function() {
    plot(cars, main = "") # here, could use main directly
    title(main = "Stopping Distance versus Speed")
}

title2 <- function() {
    plot(cars, main = "")
    title(main = list("Stopping Distance versus Speed", cex = 1.5,
              col = "red", font = 3))
}

title3 <- function() {
    ## Specifying "..." :
    plot(1, col.axis = "sky blue", col.lab = "thistle")
    title("Main Title", sub = "sub title",
          cex.main = 2,   font.main= 4, col.main= "blue",
          cex.sub = 0.75, font.sub = 3, col.sub = "red")
}

title4 <- function() {
    x <- seq(-4, 4, len = 101)
    y <- cbind(sin(x), cos(x))
    matplot(x, y, type = "l", xaxt = "n",
            main = expression(paste(plain(sin) * phi, "  and  ",
                plain(cos) * phi)),
            ylab = expression("sin" * phi, "cos" * phi), # only 1st is taken
            xlab = expression(paste("Phase Angle ", phi)),
            col.main = "blue")
    axis(1, at = c(-pi, -pi/2, 0, pi/2, pi),
         labels = expression(-pi, -pi/2, 0, pi/2, pi))
    abline(h = 0, v = pi/2 * c(-1,1), lty = 2, lwd = .1, col = "gray70")
}

plotdiff(expression(title1()), "title-1")
plotdiff(expression(title2()), "title-2")
plotdiff(expression(title3()), "title-3")
plotdiff(expression(title4()), "title-4")

plotdiffResult()
