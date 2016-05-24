##
##  i t e r s o l v e . R  Iterative Solutions of LSEs
##


itersolve <- function(A, b, x0 = NULL, 
                      nmax = 1000, tol = .Machine$double.eps^(0.5),
                      method = c("Gauss-Seidel", "Jacobi", "Richardson")) {
    stopifnot(is.numeric(A), is.numeric(b))

    n <- nrow(A)
    if (ncol(A) != n)
        stop("Argument 'A' must be a square, positive definite matrix.")
    b <- c(b)
    if (length(b) != n)
        stop("Argument 'b' must have the length 'n = ncol(A) = nrow(A).")
    if (is.null(x0)) {
        x0 <- rep(0, n)
    } else {
        stopifnot(is.numeric(x0))
        x0 <- c(x0)
        if (length(x0) != n)
            stop("Argument 'x0' must have the length 'n=ncol(A)=nrow(A).")
    }

    method <- match.arg(method)

    if (method == "Jacobi") {
        L <- diag(diag(A))
        U <- eye(n)
        beta <- 1; alpha <- 1
    } else if (method == "Gauss-Seidel") {
        L <- tril(A)
        U <- eye(n)
        beta <- 1; alpha <- 1
    } else {  # method = "Richardson"
        L <- eye(n)
        U <- L
        beta <- 0
    }

    b <- as.matrix(b)
    x <- x0 <- as.matrix(x0)
    r <- b - A %*% x0
    r0 <- err <- norm(r, "f")

    iter <- 0
    while (err > tol && iter < nmax) {
        iter <- iter + 1
        z <- qr.solve(L, r)
        z <- qr.solve(U, z)
        if (beta == 0) alpha <- drop(t(z) %*% r/(t(z) %*% A %*% z))
        x <- x + alpha * z
        r <- b - A %*% x
        err <- norm(r, "f") / r0
    }

    return(list(x = c(x), iter = iter, method = method))
}
