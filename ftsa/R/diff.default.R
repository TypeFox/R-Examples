diff.default = function (x, lag = 1, differences = 1, ...) 
{
    ismat <- is.matrix(x)
    xlen <- if (ismat) 
        dim(x)[1L]
    else length(x)
    if (length(lag) > 1L || length(differences) > 1L || lag < 
        1L || differences < 1L) 
        stop("'lag' and 'differences' must be integers >= 1")
    if (lag * differences >= xlen) 
        return(x[0])
    r <- unclass(x)
    i1 <- -1L:-lag
    if (ismat) 
        for (i in 1L:differences) r <- r[i1, , drop = FALSE] - 
            r[-nrow(r):-(nrow(r) - lag + 1), , drop = FALSE]
    else for (i in 1L:differences) r <- r[i1] - r[-length(r):-(length(r) - 
        lag + 1)]
    class(r) <- oldClass(x)
    r
}
