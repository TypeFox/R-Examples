rapca = function (x, FUN = Qn, order = 4, mean = TRUE)
{
    if (order < 1)
        stop("Order must be positive")
    X <- t(x)
    n <- nrow(X)
    p <- ncol(X)
    if (mean) {
        med <- colMeans(X)
        xx <- sweep(X, 2, med)
    }
    else xx <- X
    tmp <- La.svd(xx)
    r = sum(tmp$d > (max(n, p) * max(tmp$d) * 1e-12))
    P <- t(tmp$vt)[, 1:r]
    tmp2 <- rstep(t(xx %*% P), order = order, r = tmp$r, mean = mean)
    tmp <- P %*% tmp2$basis
    if (mean) {
        med <- c(med + tmp[, 1])
        xx <- sweep(X, 2, med)
        basis <- cbind(med, tmp[, (1:order) + 1])
        coef <- cbind(rep(1, n), xx %*% basis[, -1])
    }
    else {
        basis <- tmp
        coef <- xx %*% basis
    }
    return(list(basis = basis, coeff = coef, X = xx))
}
