"dotPlot" <-
function (x, y = 0, xlim = range(x, na.rm = TRUE), xlab = NULL, 
    scatter = FALSE, hmax = 1, base = TRUE, axes = TRUE, frame = FALSE, 
    pch = 21, pch.size = "x", labels = NULL, hcex = 1, cex = par("cex"), 
    cex.axis = par("cex.axis"), ...) 
{
    if (is.null(xlab)) 
        xlab <- deparse(substitute(x))
    x <- x[!is.na(x)]
    xpd <- par("xpd")
    par(xpd = TRUE)
    on.exit(par(xpd = xpd))
    plot(c(0, 1), c(0, 1), xlim = xlim, type = "n", axes = FALSE, 
        cex = cex, cex.axis = cex.axis, frame = frame, xlab = xlab, 
        ylab = "", ...)
    if (axes) 
        axis(1, cex.axis = cex.axis)
    if (scatter) {
        dots(x, y = y, xlim = xlim, stacked = FALSE, hmax = hmax, 
            base = base, axes = FALSE, pch = pch, pch.size = pch.size, 
            labels = labels, hcex = hcex, cex = cex, cex.axis = cex.axis)
        y <- y + 2 * strheight(pch.size, units = "user")
        xlab <- ""
        axes <- FALSE
        base = FALSE
    }
    coord <- dots(x, y, xlim = xlim, stacked = TRUE, hmax = hmax, 
        base = base, axes = FALSE, pch = pch, pch.size = pch.size, 
        labels = labels, hcex = hcex, cex = cex, cex.axis = cex.axis)
    invisible(coord)
}
