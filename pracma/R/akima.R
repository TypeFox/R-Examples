##
##  a k i m a . R  Univariate Akima Interpolation
##


akimaInterp <- function(x, y, xi) {
    stopifnot(is.numeric(x), is.numeric(y), is.numeric(xi),
               is.vector(x),  is.vector(y),  is.vector(xi))

    n <- length(x)
    if (length(y) != n)
        stop("Vectors 'x' and 'y' must be of the same length.")
    dx <- diff(x)
    if (any(dx <= 0))
        stop("Argument 'x' must be an in strictly ascending order.")
    if (any(xi < x[1]) || any(xi > x[n]))
        stop("All points in 'xi' must lie between x[1] and x[n].")

    m <- diff(y) / dx
    mm <- 2*m[1]-m[2];     mmm <- 2*mm-m[1]     # augment at left
    mp <- 2*m[n-1]-m[n-2]; mpp <- 2*mp-m[n-1]   # augment at right
    m1 <- c(mmm, mm, m, mp, mpp)                # slopes

    dm <- abs(diff(m1))
    f1 <- dm[3:(n+2)]; f2 <- dm[1:n]; f12 <- f1 + f2
    id <- which(f12 > 1e-8 * max(f12))
    b <- m1[2:n+1]
    b[id] <- (f1[id] * m1[id+1] + f2[id] * m1[id+2]) / f12[id]
    e <- (3*m - 2*b[1:n-1] - b[2:n]) / dx
    d <- (b[1:n-1] + b[2:n] - 2*m) / dx^2

    bin <- findInterval(xi,x)
    bin <- pmin(bin,n-1)
    bb <- bin[1:length(xi)]
    wj <- xi - x[bb]
    yi <- ((wj * d[bb] + e[bb]) * wj + b[bb]) * wj + y[bb]

    return(yi)
}
