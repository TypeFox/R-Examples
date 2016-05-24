###
### p o l y a r e a . R  Calculate area and center of a polygon
###


polyarea <- function(x, y) {
    if (length(x) == 0 && length(y) == 0) return(0)
    if (!(is.numeric(x) || is.complex(x)) ||
        !(is.numeric(y) || is.complex(y)))
        stop("Arguments 'x' and 'y' must be real or complex.")
    if (is.null(dim(x))) x <- matrix(x, length(x), 1)
    if (is.null(dim(y))) y <- matrix(y, length(y), 1)
    if (any(dim(x) != dim(y)))
        stop("Matrices 'x' and 'y' must be of same size.")

    n <- nrow(x); m <- ncol(x)
    z <- numeric(m)
    for (i in 1:m) {
        xi <- x[, i]
        yi <- y[, i]
        # Gauss' formula
        p1 <- sum(xi[1:(n-1)]*yi[2:n]) + xi[n]*yi[1]
        p2 <- sum(xi[2:n]*yi[1:(n-1)]) + xi[1]*yi[n]
        z[i] <- 0.5*(p1-p2)
    }
    return(z)
}


poly_center <- function(x, y) {
    stopifnot(is.numeric(x), is.numeric(y))
    n <- length(x)
    if (length(y) != n || n <= 2)
        stop("Arguments 'x' and 'y' must be of the same length >= 3.")

    parea <- polyarea(x, y)
    if (parea == 0)
        return(c(NA, NA))

    x1 <- x[1:(n-1)]; x2 <- x[2:n]
    y1 <- y[1:(n-1)]; y2 <- y[2:n]
    xy <- x1*y2 - x2*y1

    cx <- sum((x1+x2) * xy)
    cy <- sum((y1+y2) * xy)
    return(1/parea/6 * c(cx, cy))
}


poly_length <- function(x, y) {
    stopifnot(is.numeric(x), is.numeric(y))

    X  <- cbind(x, y)
    dX <- diff(X)
    return(sum(sqrt(rowSums(dX^2))))
}


poly_crossings <- function(L1, L2) {
    stopifnot(is.numeric(L1), is.numeric(L2))
    # L1, L2 marices with two rows: rbind(x, y)
    if (!is.matrix(L1) || !is.matrix(L2) || nrow(L1) != 2 || nrow(L2) != 2)
        stop("Arguments 'L1', 'L2' must be matrices with 2 rows.")

    # Utility function
    Dd <- function(x, y) {
        x1 <- x[, 1:(ncol(x)-1)]
        x2 <- x[, 2:ncol(x)]
        y1 <- repmat(as.matrix(y), 1, ncol(x1))
        y2 <- repmat(y, nrow(x2), 1)
        (x1 - y1) * (x2 - y1)
    }

    # Preliminary stuff
    x1  <- L1[1, ];    x2  <- L2[1, ]
    y1  <- L1[2, ];    y2  <- L2[2, ]
    dx1 <- diff(x1);   dy1 <- diff(y1)
    dx2 <- diff(x2);   dy2 <- diff(y2)
    n1  <- length(x1); n2  <- length(x2)

    # Determine 'signed differences'
    S1 <- dx1 * y1[1:(n1-1)] - dy1 * x1[1:(n1-1)]
    S2 <- dx2 * y2[1:(n2-1)] - dy2 * x2[1:(n2-1)]
    X1 <- outer(dx1, y2, "*") - outer(dy1, x2, "*")
    X2 <- outer(y1, dx2, "*") - outer(x1, dy2, "*")

    C1 <- Dd(X1, S1)
    C2 <- t(Dd(t(X2), S2))

    # Segments with expected intersection
    ij <- which((C1 <= 0) & (C2 <= 0), arr.ind = TRUE)
    if (length(ij) == 0)
        return(c())

    # Prepare for output
    i <- ij[, 1]; j <- ij[, 2]
    L <- dy2[j] * dx1[i] - dy1[i] * dx2[j]
    i <- i[L != 0]; j <- j[L != 0]
    L <- L[L != 0]              # avoid divisions by 0

    # Get the common points
    P <- cbind(dx2[j] * S1[i] - dx1[i] * S2[j],
               dy2[j] * S1[i] - dy1[i] * S2[j]) / cbind(L, L)
    # TO DO: throw out equal points

    colnames(P) <- c('x', 'y')
    return(P)
}
