##
##  e i g j a c o b i . R  Jacobi Eigenvalue Method
##


eigjacobi <- function(A, tol = .Machine$double.eps^(2/3)) {
    stopifnot(is.numeric(A))
    n <- nrow(A)
    if (ncol(A) != n || any(t(A) != A))
        stop("Matrix 'A' must be a square and real symmetric matrix.")

    D <- A
    V <- eye(n)

    # calculate greatest off-diagonal element
    ind <- which.max(abs(D - diag(diag(D))))
    pq <- arrayInd(ind, dim(D))
    p <- pq[1]; q <- pq[2]
    
    while (TRUE) {
        # Zero out D[p, q] and D[q, p]
        t <- D[p,q] / (D[q,q] - D[p,p])
        d <- 1/sqrt(t^2 + 1)
        s <- d * t
        R <- matrix(c(d, -s, s, d), 2, 2)
        
        D[c(p, q), ] <- t(R) %*% D[c(p, q), ]
        D[, c(p, q)] <- D[, c(p, q)] %*% R
        V[, c(p, q)] <- V[, c(p, q)] %*% R

        ind <- which.max(abs(D - diag(diag(D))))
        pq <- arrayInd(ind, dim(D))
        p <- pq[1]; q <- pq[2]
        if (abs(D[p,q]) < tol * sqrt(sum(diag(D)^2)/n))
            break
    }

    return(list(V = V, D = diag(D)))
}
