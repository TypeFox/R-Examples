##
##  p i n v . R  Pseudoinverse (Moore-Penrose Generalized Inverse)
##


pinv <- function (A, tol = .Machine$double.eps^(2/3)) {
    stopifnot(is.numeric(A), length(dim(A)) == 2, is.matrix(A))

    s <- svd(A)
    # D <- diag(s$d); Dinv <- diag(1/s$d)
    # U <- s$u; V <- s$v
    # A = U D V'
    # X = V Dinv U'

    p <- ( s$d > max(tol * s$d[1], 0) )
    if (all(p)) {
        mp <- s$v %*% (1/s$d * t(s$u))
    } else if (any(p)) {
        mp <- s$v[, p, drop=FALSE] %*% (1/s$d[p] * t(s$u[, p, drop=FALSE]))
    } else {
        mp <- matrix(0, nrow=ncol(A), ncol=nrow(A))
    }

    return(mp)
}
