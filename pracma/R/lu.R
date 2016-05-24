##
##  l u. R  LU Decomposition
##


lu <- function(A, scheme = c("kji", "jki", "ijk")) {
    stopifnot(is.numeric(A), is.matrix(A))
    n <- nrow(A)
    if (ncol(A) != n || n <= 1)
        stop("Argument 'A' must be a square, positive definite matrix.")

    scheme <- match.arg(scheme)

    if (scheme == "kji") {
        for (k in 1:(n-1)) {
            if (A[k, k] == 0)
                stop("All diagonal elements of matrix 'A' must be non-zero.")
            for (i in (k+1):n) {
                A[i, k] <- A[i, k] / A[k, k]
                A[i, (k+1):n] <- A[i, (k+1):n] - A[i, k] * A[k, (k+1):n]
            }
        }

    } else if (scheme == "jki") {
        if (A[1, 1] == 0)
            stop("All diagonal elements of matrix 'A' must be non-zero.")
        i <- 2:n
        A[i, 1] <- A[i, 1] / A[1, 1]
        for (j in 2:n) {
            if (A[j, j] == 0)
                stop("All diagonal elements of matrix 'A' must be non-zero.")
            for (k in 1:(j-1)) {
                i <- (k+1):n
                A[i, j] <- A[i, j] - A[i, k] * A[k, j]
            }
            if (j < n) {
                i <- (j+1):n
                A[i, j] <- A[i, j] / A[j, j]
            }
        }

    } else if (scheme == "ijk") {       # 'compact' Doolittle scheme
        for (i in 2:n) {                # j in 1:n
            for (j in 2:i) {
                if (A[j, j] == 0)
                    stop("All diagonal elements of matrix 'A' must be non-zero.")
                A[i, j-1] <- A[i, j-1] / A[j-1, j-1]
                k <- 1:(j-1)
                A[i, j] <- A[i, j] - A[i, k] %*% A[k, j]
            }
            if (i < n) {
                k <- 1:(i-1)
                for (j in (i+1):n) {
                    A[i, j] <- A[i, j] - A[i, k] %*% A[k, j]
                }
            }
        }
    }

    L <- eye(n) + tril(A, -1)
    U <- triu(A)
    return(list(L = L, U = U))
}


lufact <- function(A) {
    stopifnot(is.numeric(A), is.matrix(A))
    m <- nrow(A); n <- ncol(A)
    if (m != n || m == 1)
        stop("Matrix 'A' must be a square matrix with 2 rows at least.")

    detA <- 1
    rows <- 1:n
    for (p in 1:(n-1)) {
        prow <- which.max(abs(A[p:n, p])) + (p-1)
        if (p < prow) {
            rows[c(p, prow)] <- rows[c(prow, p)]
            detA <- -detA
        }
        detA <- detA * A[rows[p], p]
        if (detA == 0) {
            warning("Matrix 'A' is singular, no results computed.")
            return(list(LU = A, perm = rows, det = NA, x = NULL))
        }

        for (k in (p+1):n) {
            f <- A[rows[k], p] / A[rows[p], p]
            A[rows[k], p] <- f
            A[rows[k], (p+1):n] <- A[rows[k],(p+1):n] - f*A[rows[p], (p+1):n]
        }
    }
    detA <- detA * A[rows[n], n]
    return(list(LU = A, perm = rows, det = detA))
 }


lusys <- function(A, b) {
    stopifnot(is.numeric(A), is.matrix(A), is.numeric(b))
    b <- as.matrix(b)

    m <- nrow(A); n <- ncol(A)
    if (m != n || m == 1)
        stop("Matrix 'A' must be a square matrix with 2 rows at least.")

    x <- zeros(n, 1); y <- zeros(n, 1)
    r <- c(1:n)

    for (p in 1:(n-1)) {
        # find the pivot row for column p
        q <- which.max(abs(A[p:n, p])) + (p-1)

        # interchange rows p and q
        A[c(p, q), ] <- A[c(q, p), ]
        r[c(p, q)] <- r[c(q, p)]

        if (A[p, p] == 0)
            stop("Matrix 'A' singular: no unique solution.")

        # calculate multiplier and place
        for (k in (p+1):n) {
            a = A[k, p] / A[p, p]
            A[k, p] <- a
            A[k, (p+1):n] <- A[k, (p+1):n] - a*A[p, (p+1):n]
        }
    }
    # solve for y
    y[1] = b[r[1]]
    for (k in 2:n)
        y[k] <- b[r[k]] - A[k,1:(k-1)] %*% y[1:(k-1)]

    # solve for x
    x[n] <- y[n] / A[n, n]
    for (k in (n-1):1)
        x[k] <- (y[k] - A[k, (k+1):n] %*% x[(k+1):n]) / A[k, k]

    return(x)
}

