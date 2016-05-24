##
##  d e v a l . R
##


deval <- function(x, y, xp, idx = NULL) {
    stopifnot(is.vector(x, mode = "numeric"), is.numeric(y),
              is.vector(xp, mode= "numeric"))
    if (is.vector(y)) y <- as.matrix(y)
    if (length(x) != nrow(y))
        stop("Length of 'x' must be equal to the number of rows in 'y'.")
    if (is.unsorted(x))
        stop("Argument vector 'x' must be sorted.")
    if (is.null(idx)) idx <- 1:ncol(y)
    if (! all(idx %in% 1:ncol(y)))
        stop("Indices 'idx' must be between 1 and no. of columns of 'y'.")

    fint <- findInterval(xp, x)
    flen <- length(fint)

    yp <- matrix(NA, nrow = flen, ncol = length(idx))

    for (i in 1:flen) {
        fi <- fint[i]
        if (fi == 0) next
        if (fi < length(x)) {
            yp[i, ] <- y[fi, idx] + 
                       (xp[i] - x[fi])/(x[fi+1] - x[fi]) * (y[fi+1, idx] - y[fi, idx])
        } else {
            if (xp[i] > x[length(x)]) {
                next
            } else {
                yp[i, ] <- y[fi, idx]
            }
        }
    }

    if (flen == 1) yp <- drop(yp)
    return(yp)
}


deeve <- function(x, y, yv = 0, idx = NULL){
    stopifnot(is.vector(x, mode = "numeric"), is.numeric(y),
              is.numeric(yv), length(yv) == 1)
    if (is.vector(y)) y <- as.matrix(y)
    if (length(x) != nrow(y))
        stop("Length of 'x' must be equal to the number of rows in 'y'.")
    if (is.null(idx)) idx <- ncol(y)
    else if (length(idx) > 1) {
        idx <- idx[1]
        warning("Several indices found; only accepting the first one.")
    }

    y <- y[, idx]
    if (yv < min(y) || yv > max(y))
        return(NA)

    # findInterval() needs nondecreasingly sorted vector
    fint <- findintervals(yv, y)
    flen <- length(fint)
    if (flen == 0) return(c())

    x0 <- numeric(flen)
    for (i in 1:flen) {
        fi <- fint[i]
        if (fi < length(y)) {
            x0[i] <- (yv - y[fi]) / (y[fi+1] - y[fi]) * (x[fi+1] - x[fi]) + x[fi]
        } else {
            x0[i] <- x[fi]
        }
    }

    return(x0)
}
