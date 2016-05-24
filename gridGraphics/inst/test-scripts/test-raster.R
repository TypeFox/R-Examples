
require(grDevices)

library(gridGraphics)

raster1 <- function() {
    ## set up the plot region:
    par(bg = "thistle")
    plot(c(100, 250), c(300, 450), type = "n", xlab = "", ylab = "")
    image <- as.raster(matrix(0:1, ncol = 5, nrow = 3))
    rasterImage(image, 100, 300, 150, 350, interpolate = FALSE)
    rasterImage(image, 100, 400, 150, 450)
    rasterImage(image, 200, 300, 200 + xinch(.5), 300 + yinch(.3),
                interpolate = FALSE)
    rasterImage(image, 200, 400, 250, 450, angle = 15, interpolate = FALSE)
}

plotdiff(expression(raster1()), "raster-1")

plotdiffResult()
