`as.mefa.array` <-
function(x, ...)
{
    if (length(dim(x)) == 2)
        return(as.mefa.default(x[[1]]))
    n <- dim(x)[3]
    tmp <- list()
    tmp[[1]] <- x[,,1]
    for (i in 2:n) {
        tmp[[i]] <- tmp[[(i-1)]] + x[,,i]
    }
    m <- mefa(tmp[[n]])
    names(tmp) <- dimnames(x)[[3]]
    m$segm <- tmp
    m <- as.mefa.default(m, ...)
    m$call <- match.call()
    m
}
