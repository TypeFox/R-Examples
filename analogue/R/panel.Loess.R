panel.Loess <- function(x, y,
                        span = 1/3,
                        degree = 1,
                        family = c("symmetric", "gaussian"),
                        evaluation = 50,
                        lwd = plot.line$lwd,
                        lty = plot.line$lty,
                        col,
                        col.line = plot.line$col,
                        type,
                        ...) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 1)
        return()
    if (!missing(col)) {
        if (missing(col.line))
            col.line <- col
    }
    plot.line <- trellis.par.get("plot.line")
    smooth <- loess.smooth(x = y[ok], y = x[ok],
                           span = span,
                           family = family,
                           degree = degree,
                           evaluation = evaluation)
    panel.lines(x = smooth$y, y = smooth$x, col = col.line, lty = lty,
                lwd = lwd, ...)
}
