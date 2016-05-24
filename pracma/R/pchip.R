##
##  p c h i p . R  Piecwise Cubic Hermitean Interpolation Polynomials
##


pchip <- function(xi, yi, x) {
    stopifnot(is.numeric(xi), is.numeric(yi), is.numeric(x))
    # xi <- c(xi); yi <- c(yi); x <- c(x)
    if (!is.sorted(xi))
        stop("Argument 'xi' must be a sorted vector of real numbers.")
    n <- length(xi);
    if (length(yi) != n)
        stop("Arguments 'xi', 'yi' must be vectors of equal length.")
    if (n <= 2)
        stop("At least three points needed for cubic interpolation.")

    # First derivatives
    h <- diff(xi)
    delta <- diff(yi) / h
    d <- .pchipslopes(h, delta)

    # Piecewise polynomial coefficients
    a <- (3*delta - 2*d[1:(n-1)] - d[2:n]) / h
    b <- (d[1:(n-1)] - 2*delta + d[2:n]) / h^2;

    # Find subinterval indices k so that xi[k] <= x < xi[k+1]
    k <- rep(1, length(x))
    for (j in 2:(n-1)) {
       k[xi[j] <= x] <- j
    }

    # Evaluate interpolant
    s <- x - xi[k]
    v <- yi[k] + s*(d[k] + s*(a[k] + s*b[k]))

    return(v)
}


.pchipslopes <- function(h, delta) {

    # Slopes at interior points
    n <- length(h) + 1
    d <- numeric(length(h))
    k <- which(sign(delta[1:(n-2)]) * sign(delta[2:(n-1)]) > 0) + 1
    w1 <- 2*h[k] + h[k-1]
    w2 <- h[k]+2*h[k-1]
    d[k] <- (w1+w2) / (w1/delta[k-1] + w2/delta[k])

    # Slopes at endpoints
    d[1] <- .pchipend(h[1], h[2], delta[1], delta[2])
    d[n] <- .pchipend(h[n-1], h[n-2], delta[n-1], delta[n-2])

    return(d)
}


.pchipend <- function(h1, h2, del1, del2) {
    # Noncentered, shape-preserving, three-point formula.
    d <- ((2*h1 + h2)*del1 - h1*del2) / (h1 + h2)
    if (sign(d) != sign(del1)) {
        d <- 0
    } else if ((sign(del1) != sign(del2)) && (abs(d) > abs(3*del1))) {
        d <- 3*del1
    }
    return(d)
}


pchipfun <- function(xi, yi) {
    stopifnot(is.numeric(xi), is.numeric(yi))
    # xi <- c(xi); yi <- c(yi)
    if (!is.sorted(xi))
        stop("Argument 'xi' must be a sorted vector of real numbers.")
    n <- length(xi);
    if (length(yi) != n)
        stop("Arguments 'xi', 'yi' must be vectors of equal length.")
    if (n <= 2)
        stop("At least three points needed for cubic interpolation.")

    function(x) pchip(xi, yi, x)
}
