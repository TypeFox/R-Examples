###
### p o l y . R  Polynom
###

Poly <- function(x) {
    if (is.null(x) || length(x) == 0) return(c(1))
    if (is.vector(x, mode="numeric") || is.vector(x, mode="complex")) {
        y <- x
    } else {
        if ((is.numeric(x) || is.complex(x))  && is.matrix(x)) {
            y <- eigen(x)$values
        } else {
            stop("Argument 'x' must be a vector or square matrix.")
        }
    }

    n <- length(y)
    p <- c(1, rep(0, n))
    for (i in 1:n) {
        p[2:(i+1)] <- p[2:(i+1)] - y[i] * p[1:i]
    }
    if (all(Im(p) == 0)) p <- Re(p)
    return(p)
}
