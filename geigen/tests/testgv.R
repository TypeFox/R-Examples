testgv <- function(A,B,z) {
    V <- z$vectors         
    if( is.null(z$alpha) ) {
        ret <- A %*% V - B %*% V %*% diag(z$values)
    } else {
        ret <- A %*% V %*% diag(z$beta) - B %*% V %*% diag(z$alpha)
    }
    tol <- 100 * sqrt(.Machine$double.eps)
    all(abs(ret)<=tol)
}
