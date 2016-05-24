fd.vec <-
function (x, f, ..., eps = 1e-05) {
    n <- length(x)
    res <- matrix(0, n, n)
    ex <- pmax(abs(x), 1)
    f0 <- f(x, ...)
    for (i in 1:n) {
        x1 <- x
        x1[i] <- x[i] + eps * ex[i]
        diff.f <- c(f(x1, ...) - f0)
        diff.x <- x1[i] - x[i]
        res[, i] <- diff.f / diff.x
    }
    0.5 * (res + t(res))
}
