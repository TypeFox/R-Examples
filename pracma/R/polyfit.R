###
### p o l y f i t . R  Polynom
###

polyfit <- function(x, y, n = 1) {
    if (!is.numeric(x) || !is.numeric(y))
        stop("Arguments x and y must be numeric.")
    if (length(x) != length(y))
        stop("Vectors/matrices x and y must be of same length.")
    if (is.null(n) || n < 0 || ceiling(n) != floor(n))
        stop("Degree n must be a non-negative integer.")

    x <- x[1:length(x)]; y <- y[1:length(y)]
    A <- outer(x, seq(n, 0), "^")
    p <- qr.solve(A, y)
    return(p)
}

polyfit2 <- function(x, y, n = 1, p0 = NULL) {
    if (!is.numeric(x) || !is.numeric(y))
        stop("Arguments 'x' and 'y' must be numeric.")
    if (length(x) != length(y))
        stop("Vectors/matrices 'x' and 'y' must be of same length.")
    if (is.null(n) || n <= 0 || ceiling(n) != floor(n))
        stop("Argument 'n', order of fit, must be a positive integer.")
    if (is.null(p0))
        return(polyfit(x, y, n = n))
    else if (!is.numeric(p0) || length(p0) != 2)
        stop("Argument 'p0' must be a numeric vector of length 2.")

    x0 <- p0[1];  y0 <- p0[2]
    xx <- x - x0; yy <- y - y0

    M <- matrix(0, length(x), n)
    M[, n] <- xx
    for (i in (n-1):1) {
        M[, i] <- xx * M[, i+1]
    }
    pt <- qr.solve(M, yy)
    pt <- c(pt, y0)
    p <- numeric(n+1)
    for (i in (n+1):1) {
        p[i] <- polyval(pt, -x0)
        pt   <- polyder(pt)/(n-i+2)
    }
    return(p)
}
