`coinertiaI` <- function(X, Y, fast = TRUE) {
    A <- t(X) %*% Y
    retval <- if (fast) {
        svdA <- La.svd(A, nu = 0)
        Psi <- Y %*% t(svdA$vt)
    } else {
        svdA <- La.svd(A)
        Ksi <- X %*% svdA$u
        Psi <- Y %*% t(svdA$vt)
        L <- diag(svdA$d)^2
        retval <- list(weights = list(X = svdA$u, Y = t(svdA$vt)),
                       scores = list(X = Ksi, Y = Psi),
                       lambda = L, call = match.call())
        class(retval) <- c("coinertiaI", "fitCoinertia")
    }
    retval
}

