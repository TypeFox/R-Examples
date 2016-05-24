set_outer <-
function(X, Y, FUN = "*", ..., SIMPLIFY = TRUE, quote = FALSE)
{
    ## convenience
    nx <- deparse(substitute(X))
    if(missing(Y)) {
        Y <- X
        ny <- nx
    } else if(is.function(Y) || is.character(Y)) {
        FUN <- Y
        Y <- X
        ny <- nx
    } else ny <- deparse(substitute(Y))

    FUN <- match.fun(FUN)

    ## convert to lists
    xlabs <- LABELS(X, quote = quote)
    ylabs <- LABELS(Y, quote = quote)
    X <- as.list(X)
    Y <- as.list(Y)

    ## loop
    xrep <- rep(X, times = (ylen <- length(Y)))
    yrep <- rep(Y, each = (xlen <- length(X)))
    ret <- mapply(FUN, xrep, yrep, MoreArgs = list(...), SIMPLIFY = FALSE)

    ## simplify if sensible
    if(SIMPLIFY && all(sapply(ret, is.atomic)))
        ret <- unlist(ret, recursive = FALSE)

    ## make matrix
    dim(ret) <- c(xlen, ylen)
    dimnames(ret) <- list(xlabs, ylabs)
    ret
}

cset_outer <-
gset_outer <-
function(X, Y, FUN = "*", ..., SIMPLIFY = TRUE, quote = FALSE)
    set_outer(X = X, Y = Y, FUN = FUN, ..., SIMPLIFY = SIMPLIFY, quote = quote)
