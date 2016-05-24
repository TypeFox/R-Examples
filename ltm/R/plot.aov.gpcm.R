plot.aov.gpcm <-
function (x, line = TRUE, xlab, ylab, col, ...) {
    if (!inherits(x, "aov.gpcm") || !x$simulate.p.value)
        stop("Use only with 'aov.gpcm' objects with simulate.p.value = TRUE.\n")
    vals <- x$LRTvals
    ord <- order(vals)
    ord.vals <- vals[ord]
    df <- x$df
    z <- qchisq(ppoints(length(ord.vals)), df)
    if (missing(xlab))
        xlab <- substitute(paste("Theoretical Quantiles  ", chi[x]^2), list(x = df))
    if (missing(ylab))
        ylab <- "Sample LRT Quantiles"
    if (missing(col))
        col <- "red"
    plot(z, ord.vals, xlab = xlab, ylab = ylab, col = col, ...)
    if (line) {
        abline(a = 0, b = 1, ...)
    }
    invisible()
}
