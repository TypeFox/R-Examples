
require(grDevices) # for colours

library(gridGraphics)

image1 <- function() {
    x <- y <- seq(-4*pi, 4*pi, len = 27)
    r <- sqrt(outer(x^2, y^2, "+"))
    z <- z <- cos(r^2)*exp(-r/6)
    image(z, col  = gray((0:32)/32))
}

image2 <- function() {
    x <- y <- seq(-4*pi, 4*pi, len = 27)
    r <- sqrt(outer(x^2, y^2, "+"))
    z <- z <- cos(r^2)*exp(-r/6)
    image(z, axes = FALSE, main = "Math can be beautiful ...",
          xlab = expression(cos(r^2) * e^{-r/6}))
    contour(z, add = TRUE, drawlabels = FALSE)
}

image3 <- function() {
    # Volcano data visualized as matrix. Need to transpose and flip
    # matrix horizontally.
    image(t(volcano)[ncol(volcano):1,])
}

image4 <- function() {
    # A prettier display of the volcano
    x <- 10*(1:nrow(volcano))
    y <- 10*(1:ncol(volcano))
    image(x, y, volcano, col = terrain.colors(100), axes = FALSE)
    contour(x, y, volcano, levels = seq(90, 200, by = 5), drawlabels = FALSE,
            add = TRUE, col = "peru")
    axis(1, at = seq(100, 800, by = 100))
    axis(2, at = seq(100, 600, by = 100))
    box()
    title(main = "Maunga Whau Volcano", font.main = 4)
}

plotdiff(expression(image1()), "image-1", antialias=FALSE)
plotdiff(expression(image2()), "image-2", antialias=FALSE)
plotdiff(expression(image3()), "image-3", antialias=FALSE)
plotdiff(expression(image4()), "image-4", antialias=FALSE)

plotdiffResult()
