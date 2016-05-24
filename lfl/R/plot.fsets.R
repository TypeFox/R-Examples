plot.fsets <- function(x, ...) {
    if (!is.fsets(x)) {
        stop("'x' must be an instance of the 'fsets' class")
    }

    n <- nrow(x)
    args <- list(...)
    x <- as.list(as.data.frame(x))
    x <- llply(x, ts, start=0, frequency=n)
    x <- c(x, args)
    if (is.null(args[['xlab']])) {
        x[['xlab']] <- ''
    }
    do.call('ts.plot', x)
}
