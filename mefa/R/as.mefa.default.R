`as.mefa.default` <-
function(x, samp, taxa, ...)
{
    if (inherits(x, "mefa")) {
        ss <- if (dim(x)[3] > 1)
            melt.mefa(x) else x$xtab
        if (missing(samp))
            samp <- x$samp
        if (missing(taxa))
            taxa <- x$taxa
        y <- mefa(ss, samp, taxa, ...)
        y$call <- match.call()
        } else {
            y <- mefa(x, ...)
            y$call <- match.call()
        }
    return(y)
}

