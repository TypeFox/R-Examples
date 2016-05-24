
library(gridGraphics)

require(grDevices)

pie1 <- function() {
    pie(rep(1, 24), col = rainbow(24), radius = 0.9)
}

pie.sales <- c(0.12, 0.3, 0.26, 0.16, 0.04, 0.12)
names(pie.sales) <- c("Blueberry", "Cherry",
                      "Apple", "Boston Cream", "Other", "Vanilla Cream")

pie2 <- function() {
    pie(pie.sales) # default colours
}

pie3 <- function() {
    pie(pie.sales, col = c("purple", "violetred1", "green3",
                       "cornsilk", "cyan", "white"))
}

pie4 <- function() {
    pie(pie.sales, col = gray(seq(0.4, 1.0, length = 6)))
}

pie5 <- function() {
    pie(pie.sales, density = 10, angle = 15 + 10 * 1:6)
}

pie6 <- function() {
    pie(pie.sales, clockwise = TRUE, main = "pie(*, clockwise = TRUE)")
    segments(0, 0, 0, 1, col = "red", lwd = 2)
    text(0, 1, "init.angle = 90", col = "red")
}

pie7 <- function() {
    n <- 200
    pie(rep(1, n), labels = "", col = rainbow(n), border = NA,
        main = "pie(*, labels=\"\", col=rainbow(n), border=NA,..")
}

plotdiff(expression(pie1()), "pie-1")
plotdiff(expression(pie2()), "pie-2")
plotdiff(expression(pie3()), "pie-3")
plotdiff(expression(pie4()), "pie-4")
plotdiff(expression(pie5()), "pie-5")
plotdiff(expression(pie6()), "pie-6")
plotdiff(expression(pie7()), "pie-7")

plotdiffResult()
