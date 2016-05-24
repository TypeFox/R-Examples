fd.vec <-
function (x, f, ..., eps = 1e-06) {
    n <- length(x)
    f0 <- f(x, ...)
    res <- matrix(0, n, n)
    ex <- pmax(abs(x), 1)
    for (i in 1:n) {
        x. <- x
        x.[i] <- x[i] + eps * ex[i]
        diff.f <- c(f(x., ...) - f0)
        diff.x <- x.[i] - x[i]
        res[, i] <- diff.f / diff.x
    }
    res
}
