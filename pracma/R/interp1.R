interp1 <- function (x, y, xi = x,
             method = c("linear", "constant", "nearest", "spline", "cubic"))
{
    if (!is.vector(x, mode="numeric") || !is.vector(y, mode="numeric"))
        stop("Arguments 'x' and 'y' must be numeric vectors.")
    nx <- length(x)
    if (length(y) != nx)
        stop("Arguments 'x' and 'y' must be vectors of the same length.")
    if (nx <= 1)
        stop("Arguments 'x', 'y' must have at least a length >= 2.")

    if (min(xi) < min(x) || max(xi) > max(x))
        stop("Points 'xi' outside of range of argument 'x'.")

    e <- try(method <- match.arg(method), silent = TRUE)
    if (class(e) == "try-error") {
        warning("Unknown method: will use 'linear' interpolation.")
        method <- "linear"
    }

    if (is.unsorted(x)) {  # necessary for method 'nearest'
        warning("Points in argument in 'x' unsorted; will be sorted.")
        o <- order(x)
        x <- x[o]; y <- y[o]
    }

    if (any(duplicated(x)))
        warning("There are duplicated values in 'x'; mean will be tried.")

    if (method == "constant" || method == "linear") {
        yi <- approx(x, y, xi, method = method)$y
    } else if (method == "nearest") {
        n <- length(x)
        xx <- c(x[1], (x[2:n] + x[1:(n-1)])/2, x[n])
        yy <- c(y, y[n])
        yi <- approx(xx, yy, xi, method = "constant")$y
    } else if (method == "spline") {
        #spfun <- splinefun(x, y, method = "fmm"); yi <- spfun(xi)
        yi <- .ml.spline(x, y, xi)
    } else if (method == "cubic") {
        yi<- pchip(x, y, xi)
    } else
        stop(paste("Method", method, "not yet implemented."))

    return(yi)
}

#-- Moler's spline function --------------------------------
.ml.spline <- function(x, y, xi) {
    x <- c(x); y <- c(y)
    n <- length (x)
    if (length(y) != n) stop("Arguments 'x', 'y' must have the same length.")
    if (n < 3)          stop("spline routine: requires at least 3 points")

    # First derivatives
    h <- diff(x)
    delta <- diff(y)/h
    d <- .ml.splineslopes(h, delta)

    # Piecewise polynomial coefficients
    cc <- (3*delta - 2*d[1:(n-1)] - d[2:n])/h
    b  <- (d[1:(n-1)] - 2*delta + d[2:n])/h^2

    # Find subinterval indices k so that x(k) <= xi < x(k+1)
    m <- length(xi)
    k <- rep(1, m)
    for (j in 2:(n-1))
        k[x[j] <= xi] <- j

    # Evaluate spline
    s <- xi - x[k]
    v <- y[k] + s*(d[k] + s*(cc[k] + s*b[k]))

    return(v)
}

.ml.splineslopes <- function(h, delta) {

    # Diagonals of tridiagonal system
    n <- length(h)+1
    a <- numeric(length(h))
    b <- a; cc <- a; r <- a
    a[1:(n-2)] <- h[2:(n-1)]
    a[n-1] <- h[n-2] + h[n-1]
    b[1] <- h[2]
    b[2:(n-1)] <- 2*(h[2:(n-1)] + h[1:(n-2)])
    b[n] <- h[n-2]
    cc[1] <- h[1] + h[2]
    cc[2:(n-1)] <- h[1:(n-2)]

    # Right-hand side
    r[1] <- ((h[1] + 2*cc[1]) * h[2]*delta[1] + 
              h[1]^2 * delta[2])/cc[1]
    r[2:(n-1)] <- 3 * (h[2:(n-1)] * delta[1:(n-2)] + 
                  h[1:(n-2)] * delta[2:(n-1)])
    r[n] <- (h[n-1]^2 * delta[n-2] + 
             (2*a[n-1] + h[n-1]) * h[n-2] * delta[n-1])/a[n-1]

    # Solve tridiagonal linear system
    d <- .ml.tridisolve(a, b, cc, r)
    return(d)
}

.ml.tridisolve <- function(a, b, cc, d) {
    x <- d
    n <- length(x)

    for (j in 1:(n-1)) {
       mu <- a[j]/b[j]
       b[j+1] <- b[j+1] - mu*cc[j]
       x[j+1] <- x[j+1] - mu*x[j]
    }

    x[n] <- x[n]/b[n]
    for (j in (n-1):1)
        x[j] <- (x[j]-cc[j]*x[j+1])/b[j]

    return(x)
}
