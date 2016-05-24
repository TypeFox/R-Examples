smacofDeleteOne <- function (delta, ndim, type) {
    n <- nrow (delta)
    x <- array (0, c (n, ndim, n))
    for (i in 1:n) {
        xi <- smacofSym(delta[-i, -i], ndim = ndim, type = type)$conf
        x[((1 : n)[-i]), (1 : ndim), i] <- xi
        x[i, (1 : ndim), i] <- 0
        }
    return (x)
}