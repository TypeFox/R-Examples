require(grDevices) # for colours

library(gridGraphics)

contour1 <- function() {
    x <- -6:16
    par(mfrow = c(2, 2))
    contour(outer(x, x), drawlabels = FALSE)
    z <- outer(x, sqrt(abs(x)), FUN = "/")
    image(x, x, z)
    contour(x, x, z, col = "pink", add = TRUE, method = "edge",
            drawlabels = FALSE)
    contour(x, x, z, ylim = c(1, 6), method = "simple", 
            xlab = quote(x[1]), ylab = quote(x[2]),
            drawlabels = FALSE)
    contour(x, x, z, ylim = c(-6, 6), nlev = 20, lty = 2, main = "20 levels",
            drawlabels = FALSE)
}

contour2 <- function() {
    ## Persian Rug Art:
    x <- y <- seq(-4*pi, 4*pi, len = 27)
    r <- sqrt(outer(x^2, y^2, "+"))
    par(mfrow = c(2, 2), mar = rep(0, 4))
    for(f in pi^(0:3))
        contour(cos(r^2)*exp(-r/f),
                drawlabels = FALSE, axes = FALSE, frame = TRUE)
}

contour3 <- function() {
    rx <- range(x <- 10*1:nrow(volcano))
    ry <- range(y <- 10*1:ncol(volcano))
    ry <- ry + c(-1, 1) * (diff(rx) - diff(ry))/2
    tcol <- terrain.colors(12)
    par(pty = "s", bg = "lightcyan")
    plot(x = 0, y = 0, type = "n", xlim = rx, ylim = ry, xlab = "", ylab = "")
    u <- par("usr")
    rect(u[1], u[3], u[2], u[4], col = tcol[8], border = "red")
    contour(x, y, volcano, col = tcol[2], lty = "solid", add = TRUE,
            drawlabels = FALSE)
    title("A Topographic Map of Maunga Whau", font = 4)
    abline(h = 200*0:4, v = 200*0:4, col = "lightgray", lty = 2, lwd = 0.1)
}

# Disable antialiasing because of image() within contour1()
plotdiff(expression(contour1()), "contour-1", antialias=FALSE)
plotdiff(expression(contour2()), "contour-2")
plotdiff(expression(contour3()), "contour-3")

plotdiffResult()
