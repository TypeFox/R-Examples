##
##  t r a p z . R  Numerical integration by trapezoidal rule
##


trapz <- function(x, y) {
    if (missing(y)) {
        if (length(x) == 0) return(0)
        y <- x
        x <- seq(along=x)
    }
    if (length(x) == 0 && length(y) == 0) return(0)
    if (!(is.numeric(x) || is.complex(x)) ||
            !(is.numeric(y) || is.complex(y)) )
        stop("Arguments 'x' and 'y' must be real or complex vectors.")
    m <- length(x)
    if (length(y) != m)
        stop("Arguments 'x', 'y' must be vectors of the same length.")
    if (m <= 1) return(0.0)

    # z <- sum((x[2:m] - x[1:(m-1)]) * (y[1:(m-1)] + y[2:m]))
    # return(0.5 * z)

    xp <- c(x, x[m:1])
    yp <- c(numeric(m), y[m:1])
    n <- 2*m
    p1 <- sum(xp[1:(n-1)]*yp[2:n]) + xp[n]*yp[1]
    p2 <- sum(xp[2:n]*yp[1:(n-1)]) + xp[1]*yp[n]

    return(0.5*(p1-p2))
}


cumtrapz <- function(x, y) {
    if (missing(y)) {
        if (length(x) == 0) return(0)
        y <- x
        x <- 1:length(x)
    }
    if (length(x) == 0) return(0)
    if (!(is.numeric(x) || is.complex(x)) ||
        !(is.numeric(y) || is.complex(y)))
        stop("Arguments 'x' and 'y' must be real or complex.")

    x <- as.matrix(c(x))
    m <- length(x)
    if (is.vector(y)) y <- as.matrix(y)
    if (nrow(y) != m)
        stop("Arguments 'x' and 'y' are not compatible: nrow(y) != length(x).")

    n  <- ncol(y)
    dt <- repmat(diff(x)/2, 1, n)
    ct <- apply(dt * (y[1:(m-1), ] + y[2:m, ]), 2, cumsum)

    return(rbind(zeros(1, n), ct))
}
