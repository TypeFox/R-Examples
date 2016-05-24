##' Finds the unique commutation matrix K:
##' \eqn{K vec(A) = vec(A^t)}
##'
##' @title Finds the unique commutation matrix
##' @param m rows
##' @param n columns
##' @author Klaus K. Holst
##' @export
commutation <- function(m, n=m) {
    if (inherits(m,"matrix")) {
        n <- ncol(m)
        m <- nrow(m)
    }
    H <- function(i,j) { ## mxn-matrix with 1 at (i,j)
        Hij <- matrix(0, nrow=m, ncol=n)
        Hij[i,j] <- 1
        Hij
    }
    K <- matrix(0,m*n,m*n)
    for (i in seq_len(m))
    for (j in seq_len(n))
        K <- K + H(i,j)%x%t(H(i,j))
    K
}
