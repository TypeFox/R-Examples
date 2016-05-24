`plot.sppResponse` <- function(x, which = seq_along(x),
                               display = c("observed","fitted"),
                               xlab = "Gradient",
                               ylab = "Abundance",
                               main = NULL,
                               lcol = "red",
                               lwd = 2,
                               ...) {
    display <- match.arg(display)

    noMain <- is.null(main)
    if (noMain) {
        nams <- names(x)
    }

    ## process which - this could be logical
    if (is.logical(which))
        which <- which(which) ## yeah, really!

    for (i in which) {
        ox <- x[[i]]$observed$gradient
        oy <- x[[i]]$observed$response
        fx <- x[[i]]$fitted.values$gradient
        fy <- x[[i]]$fitted.values$response
        xlim <- range(ox, fx)
        ylim <- range(oy, fy)
        if (noMain)
            main <- nams[i]
        plot(ox, oy, xlim = xlim, ylim = ylim, ylab = ylab, xlab = xlab,
             main = main, ...)
        lines(fx, fy, col = lcol, lwd = lwd, ...)
    }

    invisible()
}
