##
##  d i s t m a t . R  Distance matrix and Hausdorff dimension
##

# If a is m x r and b is n x r then 
#
#   apply(outer(a,t(b),"-"),c(1,4),function(x)sqrt(sum(diag(x*x))))
#
# is the m x n matrix of distances between the m rows of a and 
# n rows of b.
#
# Modify, as necessary, if you want distances other than euclidean.
#
# The following code is 10-100 times faster.


distmat <- function(X, Y)
#   Computes Euclidean distance between two vectors as:
#       ||A-B|| = sqrt ( ||A||^2 + ||B||^2 - 2*A.B )
#   and vectorizes to rows in two matrices (or vectors).
{
    if (!is.numeric(X) || !is.numeric(Y))
        stop("X and Y must be numeric vectors or matrices.")
    if (is.vector(X)) dim(X) <- c(1,length(X))
    if (is.vector(Y)) dim(Y) <- c(1,length(Y))
    if (ncol(X) != ncol(Y))
        stop("X and Y must have the same number of columns.")

    m  <- nrow(X); n <- nrow(Y)
    XY <- X %*% t(Y)    # (m,n)-matrix
    XX <- matrix( rep(apply(X*X, 1, sum), n), m, n, byrow=F )
    YY <- matrix( rep(apply(Y*Y, 1, sum), m), m, n, byrow=T )

    sqrt(pmax(XX + YY - 2*XY, 0))
}


pdist <- function(X) {
    distmat(X, X)
}


pdist2 <- function(X, Y) {
    distmat(X, Y)
}


hausdorff_dist <- function(P, Q) {
    stopifnot(is.numeric(P), is.numeric(Q))
    if (is.vector(P)) P <- matrix(P, ncol = 1)
    if (is.vector(Q)) Q <- matrix(Q, ncol = 1)
    if (ncol(P) != ncol(Q))
        stop("'P' and 'Q' must have the same number of columns.")

    D <- distmat(P, Q)

    # directional Hausdorff dimension
    dhd_PQ <- max(apply(D, 1, min))
    dhd_QP <- max(apply(D, 2, min))

    return(max(dhd_PQ, dhd_QP))
}
