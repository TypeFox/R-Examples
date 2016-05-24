interleave <- function (..., append.source = TRUE, sep = ": ", drop = FALSE)
{
    ## this function has been copied from package gdata,
    ## version 2.4.2
    sources <- list(...)
    sources[sapply(sources, is.null)] <- NULL
    sources <- lapply(sources, function(x) if (is.matrix(x) ||
        is.data.frame(x))
        x
    else t(x))
    nrows <- sapply(sources, nrow)
    mrows <- max(nrows)
    if (any(nrows != mrows & nrows != 1))
        stop("Arguments have differening numbers of rows.")
    sources <- lapply(sources, function(x) if (nrow(x) == 1)
        x[rep(1, mrows), , drop = drop]
    else x)
    tmp <- do.call("rbind", sources)
    nsources <- length(sources)
    indexes <- outer((0:(nsources - 1)) * mrows, 1:mrows, "+")
    
    retval <- tmp[indexes, , drop = drop]
    if (append.source && !is.null(names(sources)))
        if (!is.null(row.names(tmp)))
            row.names(retval) <- paste(format(row.names(retval)),
                format(names(sources)), sep = sep)
        else row.names(retval) <- rep(names(sources), mrows)
    retval
}
