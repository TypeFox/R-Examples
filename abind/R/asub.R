asub <- function(x, idx, dims=seq(len=max(length(dim(x)), 1)), drop=NULL) UseMethod("asub")
asub.default <- function(x, idx, dims=seq(len=max(length(dim(x)), 1)), drop=NULL) {
    # Do arbitrary indexing of x as positions in dims
    if (length(dims)>1 && !is.list(idx))
        stop("idx must be a list when length dims>1")
    if (!is.list(idx))
        idx <- list(idx)
    if (length(idx) != length(dims))
        stop("idx has different number of indices than dims")
    # Construct a skeleton call
    xic <- Quote(x[,drop=drop])
    d <- dim(x)
    if (is.null(d))
        d <- length(x)
    if (any(dims < 1 | dims > length(d)))
        stop("dims out of range")
    # Now duplicate the empty index argument the appropriate number of times
    xic <- xic[c(1, 2, rep(3, length(d)), 4)]
    if (is.null(drop)) {
        xic <- xic[-length(xic)]
    } else {
        xic[[length(xic)]] <- drop
    }
    for (i in seq(along=dims))
        if (!is.null(idx[[i]]))
            xic[2+dims[i]] <- idx[i]
    return(eval(xic)) # , envir=parent.frame(), enclos=baseenv()))
}
