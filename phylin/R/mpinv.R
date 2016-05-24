mpinv <-
function(A, eps = 1e-13) {
    #Check also MASS package functions for the inverse of a matrix
    s <- svd(A)
    e <- s$d
    e[e > eps] <- 1/e[e > eps]
    return(s$v %*% diag(e) %*% t(s$u))
}
