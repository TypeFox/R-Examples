
require(grDevices) # for colours

library(gridGraphics)

set.seed(1)
tN <- table(Ni <- stats::rpois(100, lambda = 5))

barplot1 <- function() {
    r <- barplot(tN, col = rainbow(20))
    #- type = "h" plotting *is* 'bar'plot
    lines(r, tN, type = "h", col = "red", lwd = 2)
}

barplot2 <- function() {
    barplot(tN, space = 1.5, axisnames = FALSE,
            sub = "barplot(..., space= 1.5, axisnames = FALSE)")
}

barplot3 <- function() {
    barplot(VADeaths, beside = TRUE)
}

barplot4 <- function() {
    mp <- barplot(VADeaths) # default
    tot <- colMeans(VADeaths)
    text(mp, tot + 3, format(tot), xpd = TRUE, col = "blue")
}

barplot5 <- function() {
    barplot(VADeaths, beside = TRUE,
            col = c("lightblue", "mistyrose", "lightcyan",
                "lavender", "cornsilk"),
            legend = rownames(VADeaths), ylim = c(0, 100))
    title(main = "Death Rates in Virginia", font.main = 4)
}

barplot6 <- function() {
    hh <- t(VADeaths)[, 5:1]
    mybarcol <- "gray20"
    mp <- barplot(hh, beside = TRUE,
                  col = c("lightblue", "mistyrose",
                      "lightcyan", "lavender"),
                  legend = colnames(VADeaths), ylim = c(0,100),
                  main = "Death Rates in Virginia", font.main = 4,
                  sub = "Faked upper 2*sigma error bars", col.sub = mybarcol,
                  cex.names = 1.5)
    segments(mp, hh, mp, hh + 2*sqrt(1000*hh/100), col = mybarcol, lwd = 1.5)
    mtext(side = 1, at = colMeans(mp), line = -2,
          text = paste("Mean", formatC(colMeans(hh))), col = "red")
}

barplot7 <- function() {
    # Bar shading example
    barplot(VADeaths, angle = 15+10*1:5, density = 20, col = "black",
            legend = rownames(VADeaths))
    title(main = list("Death Rates in Virginia", font = 4))
}

barplot8 <- function() {
    # border :
    barplot(VADeaths, border = "dark blue") 
}

barplot9 <- function() {
    # log scales (not much sense here):
    barplot(tN, col = heat.colors(12), log = "y")
}

barplot10 <- function() {
    barplot(tN, col = gray.colors(20), log = "xy")
}

barplot11 <- function() {
    # args.legend
    barplot(height = cbind(x = c(465, 91) / 465 * 100,
                y = c(840, 200) / 840 * 100,
                z = c(37, 17) / 37 * 100),
            beside = FALSE,
            width = c(465, 840, 37),
            col = c(1, 2),
            legend.text = c("A", "B"),
            args.legend = list(x = "topleft"))
}

plotdiff(expression(barplot1()), "barplot-1")
plotdiff(expression(barplot2()), "barplot-2")
plotdiff(expression(barplot3()), "barplot-3")
plotdiff(expression(barplot4()), "barplot-4")
plotdiff(expression(barplot5()), "barplot-5")
plotdiff(expression(barplot6()), "barplot-6")
plotdiff(expression(barplot7()), "barplot-7")
plotdiff(expression(barplot8()), "barplot-8")
plotdiff(expression(barplot9()), "barplot-9")
plotdiff(expression(barplot10()), "barplot-10", width=14)
plotdiff(expression(barplot11()), "barplot-11")

plotdiffResult()
