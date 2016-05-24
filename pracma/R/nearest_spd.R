nearest_spd <- function(A) {
    stopifnot(is.numeric(A), is.matrix(A))
    eps <- .Machine$double.eps

    m <- nrow(A); n <- ncol(A)
    if (m != n) {
        stop("Argument 'A' must be a square matrix.")
    } else if (n == 1 && A <= 0)
        return(as.matrix(eps))

    B <- (A + t(A)) / 2                 # symmetrize A
    svdB <- svd(B)                      # H is symmetric polar factor of B
    H <- svdB$v %*% diag(svdB$d) %*% t(svdB$v)

    Ahat <- (B + H) / 2
    Ahat <- (Ahat + t(Ahat)) / 2

    # Test that Ahat is in fact positive-definite;
    # if it is not so, then tweak it just a bit.
    k <- 0; not_pd <- TRUE
    while (not_pd) {
        k <- k + 1
        try_R <- try(chol(Ahat), silent = TRUE)
        if (class(try_R) == "try-error") {
            mineig <- min(eigen(Ahat, symmetric = TRUE, only.values = TRUE)$values)
            Ahat = Ahat + (-mineig*k^2 + eps(mineig)) * diag(1, n)
        } else
            not_pd <- FALSE
    }
    Ahat
}
