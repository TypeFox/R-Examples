##
##  c u b i c s p l i n e . R  Interpolating Cubic Spline
##


cubicspline <- function(x, y, xi = NULL, endp2nd = FALSE, der = c(0, 0)) {
    n <- length(x)
    h <- x[2:n] - x[1:(n-1)]
    e <- 2 * c(h[1], h[1:(n-2)] + h[2:(n-1)], h[n-1])
    A <- Diag(e) + Diag(h, -1) + Diag(h, 1)
    d <- (y[2:n] - y[1:(n-1)]) / h
    rhs <- 3* (d[2:(n-1)] - d[1:(n-2)])
    der0 <- der[1]; dern <- der[2]
    if (endp2nd) {
        A[1, 1] <- 2 * h[1];   A[1, 2] <- h[1]
        A[n, n] <- 2 * h[n-1]; A[n-1, n-2] <- h[n-1]
        rhs <- c(3*(d[1] - der0), rhs, 3*(dern - d[n-1]))
    } else {
        A[1, ] <- 0; A[1, 1] <- 1
        A[n, ] <- 0; A[n, n] <- 1
        rhs <- c(der0, rhs, dern)
    }
    S <- zeros(n, 4)
    S[, 3] <- solve(A, rhs)
    for (m in 1:(n-1)) {
        S[m,4] = (S[m+1,3]-S[m,3]) / 3 / h[m]
        S[m,2] = d[m] - h[m]/3 * (S[m + 1,3] + 2*S[m,3])
        S[m,1] = y[m]
    }
    S <- S[1:(n-1), 4:1]
    pp <- mkpp(x, S)
    
    if (is.null(xi)) {
        return(pp)
    } else {
        yi <- ppval(pp,xi)
        return(yi)
    }
}
