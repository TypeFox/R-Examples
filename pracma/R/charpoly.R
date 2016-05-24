charpoly <- function(a, info = FALSE) {
    stopifnot(is.numeric(a), is.matrix(a))
    n <- nrow(a); m <- ncol(a)
    if (n != m || n < 2)
        stop("Argument 'a' must be a square matrix.")
    if (n > 100)
        cat("The algorithm will be *very* slow for n > 100.\n")

    p <- rep(1, n+1)

    a1 <- a
    for (k in 2:n) {
        p[k] <- -1 * sum(diag(a1))/(k-1)
        if (k == n) a2 <- a1
        a1 <- a %*% (a1 + p[k] * diag(1, n))
    }
    p[n+1] <- -1 * sum(diag(a1))/n

    if (info) {
        adet <- (-1)^n * p[n+1]
        if (adet != 0)
            ainv <- -1 * (a2 + p[n] * diag(1, n))/p[n+1]
        else
            ainv = NaN

        # determine accuracy of the computation
        e <- a2 %*% a + p[n] *a - adet * diag(1, n)
        e <- max(abs(e))
        cat("Error term:", e, "\n")
    }

    if (info) return(list(cp = p, det = adet, inv = ainv))
    else      return(p)
}
