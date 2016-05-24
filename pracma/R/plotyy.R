##
##  p l o t y y . R  Two-ordinates Plot
##


plotyy <- function( x1, y1, x2, y2, gridp = TRUE, box.col = "grey",
                                    type = "l", lwd = 1, lty = 1,
                                    xlab = "x", ylab = "y", main = "",
                                    col.y1 = "navy", col.y2 = "maroon", ...)
{
    stopifnot(is.numeric(x1), is.numeric(y1), is.numeric(x2), is.numeric(y2))

    y1pretty <- pretty(y1); y1l <- min(y1pretty); y1u <- max(y1pretty)
    y2pretty <- pretty(y2); y2l <- min(y2pretty); y2u <- max(y2pretty)

    ptrans <- function(y) y1l + (y - y2l)/(y2u - y2l) * (y1u - y1l)

    y1pretty <- pretty(c(y1, ptrans(y2)))

    opar <- par(mar = c(4.1, 4.1, 3.1, 3.1))
    plot(range(c(x1, x2)), range(y1pretty),
        xlab = xlab, ylab = ylab, main = main,
        type = "n", yaxt = "n", bty = "n", ...)
    box(col = box.col)

    mx <- axis(side = 2, at = y1pretty,  labels = FALSE, col = col.y1)
    my <- axis(side = 4, at = ptrans(y2pretty), labels = FALSE, col = col.y2)
    mtext(mx, side = 2, line = 1, at = mx, col = col.y1)
    mtext(y2pretty, side = 4, line = 1, at = my, col = col.y2)


    if (gridp) grid()
    points(x1, y1,  type = type, col = col.y1, lwd = lwd, lty = lty)
    points(x2, ptrans(y2), type = type, col = col.y2, lwd = lwd, lty = lty)
    par(opar)


    invisible()
}
