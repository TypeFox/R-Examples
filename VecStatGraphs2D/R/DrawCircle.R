DrawCircle <- function (x, y, radius, nv = 100, border = NULL, col = NA, lty = 1, lwd = 1) 
{
    angle.inc <- 2 * pi/nv
    angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
    if (length(col) < length(radius)) 
        col <- rep(col, length.out = length(radius))
    for (circle in 1:length(radius)) {
        xv <- cos(angles) * radius[circle] + x
        yv <- sin(angles) * radius[circle] + y
        polygon(xv, yv, border = border, col = col[circle], lty = lty, 
            lwd = lwd)
    }
    invisible(list(x = xv, y = yv))
}
