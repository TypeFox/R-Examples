`as.mefa.list` <-
function(x, ...)
{
    n <- length(x)
    if (isS4(x[[1]]))
        for (i in 1:n) {
            x[[i]] <- Matrix::as.matrix(x[[i]])
        }
    tmp <- x
    if (n == 1 && length(dim(tmp[[1]])) == 2)
        return(as.mefa.default(tmp[[1]]))
    for (i in 2:n) {
        tmp[[i]] <- tmp[[(i-1)]] + x[[i]]
    }
    m <- mefa(tmp[[n]])
    m$segm <- tmp
    m <- as.mefa.default(m, ...)
    m$call <- match.call()
    m
}
