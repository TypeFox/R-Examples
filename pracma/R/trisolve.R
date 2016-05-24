##
##  t r i s o l v e . R  Tridiagonal Linear System Solver
##


trisolve <- function(a, b, d, rhs) {
    stopifnot(is.numeric(a), is.numeric(b), is.numeric(d), is.numeric(rhs))

    n <- length(a)
    if (n < 3)
        stop("Argument 'n' must be greater or equal 3.")
    if (length(b) != n-1 || length(d) != n-1)
        stop("Vectors 'b' and 'd' must be of a length of length(a)-1.")
    if (length(rhs) != n)
        stop("Vector 'rhs' must be of the same length as 'a'.")
    else
        x <- rhs

    b <- c(b, 0)
    for (i in 1:(n-1)) {
        if (d[i] != 0) {
            t <- a[i]/d[i]
            si <- 1/sqrt(1+t*t)
            co <- t*si
            a[i] <- a[i]*co + d[i]*si

            h <- b[i]
            b[i] <- h*co + a[i+1]*si
            a[i+1] <- -h*si + a[i+1]*co
            d[i] <- b[i+1]*si
            b[i+1] <- b[i+1]*co

            h <- x[i]
            x[i] <- h*co + x[i+1]*si
            x[i+1] <- -h*si + x[i+1]*co
        }
    }

    if (any(a == 0.0))
        stop("Triangular matrix is singular -- system not solvable.")

    x[n]   <- x[n]/a[n]
    x[n-1] <- ( x[n-1] - b[n-1]*x[n] ) / a[n-1]
    for (i in (n-2):1) {
        x[i] <- ( x[i] - b[i]*x[i+1] - d[i]*x[i+2] ) / a[i]
    }
    return(x)
}
