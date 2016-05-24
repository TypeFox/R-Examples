
require(graphics)

library(gridGraphics)

palette1 <- function() {
    palette(gray(seq(0,.9,len = 25))) # gray scales
    matplot(outer(1:100, 1:30), type = "l", lty = 1,lwd = 2, col = 1:30,
            main = "Gray Scales Palette",
            sub = "palette(gray(seq(0, .9, len=25)))")
}

palette2 <- function() {
    ## on a device where alpha-transparency is supported,
    ## use 'alpha = 0.3' transparency with the default palette :
    palette(gray(seq(0,.9,len = 25))) # gray scales
    mycols <- adjustcolor(palette(), alpha.f = 0.3)
    opal <- palette(mycols)
    set.seed(1)
    x <- rnorm(1000); xy <- cbind(x, 3*x + rnorm(1000))
    plot (xy, lwd = 2,
          main = "Alpha-Transparency Palette\n alpha = 0.3")
    xy[,1] <- -xy[,1]
    points(xy, col = 8, pch = 16, cex = 1.5)
}

palette3 <- function() {
    palette(rainbow(10))
    plot(1:10, pch=16, col=1:10, cex=3)
    palette(heat.colors(10))
    points(10:1, pch=16, col=1:10, cex=3)
}

plotdiff(expression(palette1()), "palette-1")
plotdiff(expression(palette2()), "palette-2")
plotdiff(expression(palette3()), "palette-3")

plotdiffResult()
