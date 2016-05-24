`screeplot.symcoca` <- function(x, type = "b",
                                xlab = NULL, ylab = NULL, ...) {
    if (is.null(ylab)) {
        ylab <- "Eigenvalue"
    }
    if (is.null(xlab)) {
        xlab <- "Co-CA Axis"
    }
    evals <- eigenvals(x)
    xvals <- seq_along(evals)
    plot(xvals, evals, type = type, xlab = xlab, ylab = ylab, ...)
    invisible()
}
