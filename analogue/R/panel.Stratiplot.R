`panel.Stratiplot` <- function(x, y, type = "l",
                               col,
                               pch = plot.symbol$pch,
                               cex = plot.symbol$cex,
                               col.line = plot.line$col,
                               col.symbol = plot.symbol$col,
                               col.refline = ref.line$col,
                               col.smooth = "red",
                               col.poly = plot.line$col,
                               lty = plot.line$lty,
                               lwd = plot.line$lwd,
                               lty.smooth = plot.line$lty,
                               lwd.smooth = 2,
                               lwd.h = 3,
                               fill = plot.symbol$fill,
                               zones = NULL,
                               col.zones = plot.line$col,
                               lty.zones = plot.line$lty,
                               lwd.zones = plot.line$lwd,
                               gridh = -1, gridv = -1,
                               ...) {
    if (all(is.na(x) | is.na(y)))
        return()
    x <- as.numeric(x)
    y <- as.numeric(y)
    plot.symbol <- trellis.par.get("plot.symbol")
    plot.line <- trellis.par.get("plot.line")
    ref.line <- trellis.par.get("reference.line")
    if (!missing(col)) {
        if (missing(col.line))
            col.line <- col
        if (missing(col.symbol))
            col.symbol <- col
    }
    panel.refline(v = 0, col.line = ref.line$col, ...)
    if ("o" %in% type || "b" %in% type)
        type <- c(type, "p", "l")
    if ("g" %in% type)
        panel.grid(h = gridh, v = gridv, col.line = col.refline, ...)
    if("l" %in% type)
        panel.lines(x = x, y = y, col = col.line,
                    lty = lty, lwd = lwd, ...)
    if("p" %in% type)
        panel.points(x = x, y = y, cex = cex, fill = fill,
                     col = col.symbol, pch = pch, ...)
    if ("h" %in% type) {
        panel.lines(x = x, y = y, type = "H", col = col.line, lty = lty,
                    lwd = lwd.h, lineend = "butt", ...)
    }
    if("poly" %in% type)
        panel.polygon(x = c(0, x, 0), y = c(y[1], y, y[length(y)]),
                      border = col.poly, col = col.poly, ...)
    if("smooth" %in% type)
        panel.Loess(x, y, col = col.smooth, lwd = lwd.smooth,
                    lty = lty.smooth, ...)
    if(!is.null(zones) && is.numeric(zones))
        panel.abline(h = zones, col = col.zones, lwd = lwd.zones,
                     lty = lty.zones, ...)
}
