##
##  i n t e r p 2 . R  2-D Interpolation
##


interp2 <- function(x, y, Z, xp, yp,
                    method = c("linear", "nearest", "constant")) {
    stopifnot(is.numeric(x),  is.numeric(y),  is.numeric(Z),
              is.numeric(xp), is.numeric(yp), is.matrix(Z) )

    lx <- length(x); ly <- length(y)
    if (ncol(Z) != lx || nrow(Z) != ly)
        stop("Required: 'length(x) = ncol(Z)' and 'length(y) = nrow(Z)'.")
    n  <- length(xp)
    if (length(yp) != n)
        stop("Length of vectors 'xp' and 'yp' must be the same.")

    method <- match.arg(method)

    Z <- t(Z)         # TODO: Exchange Z[x, y] with Z[y, x] instead!
    v  <- numeric(n)

    xi <- findInterval(xp, x)
    yi <- findInterval(yp, y)

    i0 <- which(xp < min(x) | xp > max(x) | yp < min(y) | yp > max(y))

    if (method == "linear") {
        for (k in 1:length(xp)) {
            if ( k %in% i0) {
                v[k] <- NA
            } else {
                i <- xi[k]; j <- yi[k]
                if (i == lx) i <- i-1
                if (j == ly) j <- j-1
                A <- matrix(c(1, x[i],   y[j],   x[i]*y[j],
                              1, x[i],   y[j+1], x[i]*y[j+1],
                              1, x[i+1], y[j],   x[i+1]*y[j],
                              1, x[i+1], y[j+1], x[i+1]*y[j+1]),
                            nrow = 4, ncol = 4, byrow = TRUE)
                b <- Z[cbind(c(i,i,i+1,i+1), c(j,j+1,j,j+1))]
                v[k] <- sum(solve(A, b) * c(1, xp[k], yp[k], xp[k]*yp[k]))
            } 
        }

    } else if (method == "constant") {
        for (k in 1:length(xp)) {
            if ( k %in% i0) {
                v[k] <- NA
            } else {
                v[k] <- Z[xi[k], yi[k]]
            }
        }

    } else if (method == "nearest") {
        for (k in 1:length(xp)) {
            if (k %in% i0) {
                v[k] <- NA
            } else {
                i <- xi[k]
                if (i == lx) i <- i-1
                j <- yi[k]
                if (j == ly) j <- j-1
                if (xp[k] <= (x[i] + x[i+1])/2) {
                    if (yp[k] <= (y[j] + y[j+1])/2) {
                        v[k] <- Z[i, j]
                    } else {
                        v[k] <- Z[i, j+1]
                    }
                } else {
                    if (yp[k] <= (y[j] + y[j+1])/2) {
                        v[k] <- Z[i+1, j]
                    } else {
                        v[k] <- Z[i+1, j+1]
                    }
                }
            }
        }

    } else
        stop("Method 'cubic' and others are not yet implemented.")
    
    return(v)
}
