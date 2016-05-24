##
##  l 1 l i n r e g . R
##


L1linreg <- function(A, b, p = 1, tol = 1e-07, maxiter = 200) {
    stopifnot(is.numeric(A), is.numeric(b))

    if (is.vector(A)) A <- as.matrix(A)
    n <- nrow(A); m <- ncol(A)
    b <- c(b)
    if (length(b) != n)
        stop("Arguments 'A' and 'b' are not compatible.")
    if (p > 1)
            stop("Argument 'p' must be smaller or equal to 1.")

    B <- Bold <- qr.solve(A, b)
    Bold[1] <- Bold[1] + 10*tol

    iter <- 1
    while (max(abs(B - Bold)) > tol && iter <= maxiter) {
        Bold <- B
        e <- c( pmax(abs(A %*% Bold - b), tol) )
        if (p == 1) {
            w <- sqrt(1/e)
        } else {
            w <- sqrt(e^(p-2))
        }
        B <- c(qr.solve(diag(w) %*% A, w * b))
        iter <- iter + 1
    }
    if (iter > maxiter)
        warning("Reached max. no. of iterations; may not have converged.")

    return(list(x = B, reltol = max(abs(B - Bold)), niter = iter))
}
