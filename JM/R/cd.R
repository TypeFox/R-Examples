cd <-
function (x, f, ..., eps = 1e-03) {
    n <- length(x)
    res <- numeric(n)
    ex <- pmax(abs(x), 1)
    for (i in 1:n) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- x1[i] - x2[i]
        res[i] <- diff.f / diff.x
    }
    res
}
