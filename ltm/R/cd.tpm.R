cd.tpm <-
function (x, f, ..., k, eps = 1e-04) {
    res <- numeric(k)
    ex <- pmax(abs(x), 1)
    for (i in 1:k) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- 2 * max(abs(c(x1[i] - x[i], x2[i] - x[i])))
        res[i] <- diff.f / diff.x
    }
    res
}
