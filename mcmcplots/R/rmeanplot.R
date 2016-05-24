rmeanplot <- function (mcmcout, parms=NULL, regex=NULL, random=NULL, leaf.marker="[\\[_]", ylim=NULL, auto.layout=TRUE, mar=c(2.0, 2.0, 1.5, 0.25) + 0.1, col=NULL, lty=1, plot.title = NULL, main=NULL, greek = FALSE, style=c("gray", "plain"), ...) {
    mcmcout <- convert.mcmc.list(mcmcout)
    if (is.null(varnames(mcmcout))) {
        warning("Argument 'mcmcout' did not have valid variable names, so names have been created for you.")
        varnames(mcmcout) <- varnames(mcmcout, allow.null = FALSE)
    }
    parnames <- parms2plot(varnames(mcmcout), parms, regex, random, leaf.marker)
    if (length(parnames) == 0)
        stop("No parameters matched arguments 'parms' or 'regex'.")
    mcmcout <- lapply(mcmcout, function(mco) apply(mco[, parnames, drop=FALSE], 2, function(x) cumsum(x)/seq_along(x)))
    traplot(mcmcout, parms=NULL, regex=NULL, random, leaf.marker, ylim, auto.layout, mar, col, lty, plot.title, main, greek, style, ...)
}

