plot.hrsize <- function (x, ...)
{
    ## Verifications
    if (!inherits(x, "hrsize"))
        stop("should be of class hrsize")

    ## Graphical settings
    opar <- par(mfrow = n2mfrow(ncol(x)))
    on.exit(par(opar))

    ## The labels
    if (!is.null(attr(x, "xlabel"))) {
       xlabel <- attr(x, "xlabel")
    } else {
       xlabel <- "Home-range level"
    }
    if (!is.null(attr(x, "ylabel"))) {
        ylabel <- attr(x, "ylabel")
    } else {
        ylabel <- "Home-range size"
    }

    ## The plot
    for (i in 1:ncol(x)) {
        plot(as.numeric(row.names(x)), x[, i], main = names(x)[i],
            pch = 16, cex = 0.5, xlab = xlabel, ylab = ylabel)
        lines(as.numeric(row.names(x)), x[, i])
    }
}
