lsqnonneg <- function(C, d) {
    stopifnot(is.numeric(C), is.numeric(d))
    if (!is.matrix(C) || !is.vector(d))
        stop("Argument 'C' must be a matrix, 'd' a vector.")
    m <- nrow(C); n <- ncol(C)
    if (m != length(d))
        stop("Arguments 'C' and 'd' have nonconformable dimensions.")

    tol = 10 * eps() * norm(C, type = "2") * (max(n, m) + 1)

    x  <- rep(0, n)             # initial point
    P  <- logical(n); Z <- !P   # non-active / active columns
    
    resid <- d - C %*% x
    w <- t(C) %*% resid
    wz <- numeric(n)

    # iteration parameters
    outeriter <- 0; it <- 0
    itmax <- 3 * n; exitflag <- 1

    while (any(Z) && any(w[Z] > tol)) {
        outeriter <- outeriter + 1
        z <- numeric(n)
        wz <- rep(-Inf, n)
        wz[Z] <- w[Z]
        im <- which.max(wz)
        P[im] <- TRUE; Z[im] <- FALSE
        z[P] <- qr.solve(C[, P], d)
        
        while (any(z[P] <= 0)) {
            it <- it + 1
            if (it > itmax) stop("Iteration count exceeded.")

            Q <- (z <= 0) & P
            alpha <- min(x[Q] / (x[Q] - z[Q]))
            x <- x + alpha*(z - x)
            Z <- ((abs(x) < tol) & P) | Z
            P <- !Z
            z <- numeric(n)
            z[P] <- qr.solve(C[, P], d)
        }
    x <- z
    resid <- d - C %*% x
    w <- t(C) %*% resid
    }
    return(list(x = x, resid.norm = sum(resid*resid)))
}
