`plot.evenSample` <- function(x, add = FALSE, xlim = NULL, ylim = NULL,
                              col = "grey", border = "grey", lty = NULL,
                              ylab, xlab, main = NULL, sub = NULL,
                              ann = TRUE, axes = TRUE,
                              ...) {
    brks <- attr(x, "breaks")
    y <- unclass(x)
    nb <- length(brks)
    if (is.null(xlim)) {
        xlim <- range(brks)
    }
    dev.hold()
    on.exit(dev.flush())
    if (!add) {
        if (is.null(ylim)) {
            ylim <- range(y, 0)
        }
        if (missing(ylab)) {
            ylab <- "Number of Samples"
        }
        if (missing(xlab)) {
            xlab <- attr(y, "gradient")
        }
        plot.new()
        plot.window(xlim, ylim, "")
        if (ann) {
            title(main = main, ylab = ylab, xlab = xlab, sub = sub, ...)
        }
        if (axes) {
            axis(side = 1, ...)
            axis(side = 2, ...)
        }
    }
    rect(brks[-nb], 0, brks[-1L], y, col = col, border = border,
         lty = lty)
    box(...)
    invisible()
}
