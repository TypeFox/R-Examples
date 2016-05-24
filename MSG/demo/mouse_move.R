xx = runif(100)
yy = runif(100)
mousemove <- function(buttons, x, y) {
    r = 0.2
    idx = (x - r < xx & xx < x + r) & (y - r < yy & yy < y +
        r)
    plot(xx, yy, type = "n")
rect(-1,-1,2,2,col='black')
# a rectangle around the cursor
    rect(x - r, y - r, x + r, y + r, border = "yellow", lty = 2)
    points(xx[idx], yy[idx], col = "yellow", cex = 2)
    points(xx[!idx], yy[!idx], col = "red")
    NULL
}
mousedown <- function(buttons, x, y) {
    "Done"
}
par(mar = rep(0, 4), pch = 20)
if (interactive()) {
    plot(xx, yy, type = "n")
    getGraphicsEvent("Click mouse to exit", onMouseDown = mousedown,
                     onMouseMove = mousemove)
}
mousemove(, 0.2, 0.3)
arrows(0.22, 0.27, 0.18, 0.33, 0.15, lwd = 4, col = "white")
