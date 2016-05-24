expand.table <-
function (x) 
{
    if (is.null(dimnames(x)) == TRUE) 
        stop("must have dimnames")
    if (any(names(dimnames(x)) == "")) 
        stop("must have names")
    tablevars <- expand.grid(rev(dimnames(x)))
    if (length(dim(x)) > 1) {
        ftablex <- ftable(x)
        counts <- as.vector(t(ftablex[, 1:ncol(ftablex)]))
    } else {
        counts <- as.vector(x)
    }
    expansion.index <- rep(1:nrow(tablevars), counts)
    newdat <- tablevars[expansion.index, , drop=FALSE]
    row.names(newdat) <- 1:nrow(newdat)
    revnames <- rev(names(newdat))
    newdat[, revnames, drop=FALSE]
}
