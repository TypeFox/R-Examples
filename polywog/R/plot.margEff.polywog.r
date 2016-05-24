##' Plot marginal effects
##'
##' Generates density plots of the observationwise marginal effects computed by
##' \code{\link{margEff.polywog}}.
##' @param x output of \code{\link{margEff.polywog}}.
##' @param ... plotting parameters to be passed to \code{\link{plot.density}}.
##' @return Data frame containing the variables whose densities were plotted,
##' invisibly.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @method plot margEff.polywog
##' @export
plot.margEff.polywog <- function(x, ...)
{
    if (!is.list(x)) {
        x <- list(x)
        names(x) <- attr(x[[1]], "xvar")
    }

    x <- MEtoDF(x)

    ## Obtain dimensions
    mfrow <- ceiling(sqrt(length(x)))
    mfcol <- ceiling(length(x) / mfrow)
    op <- par(mfrow = c(mfrow, mfcol))
    on.exit(par(op))

    ## Plot each density plot
    for (i in seq_along(x))
        plot(density(x[[i]]), main = names(x)[i], ...)

    invisible(x)
}
