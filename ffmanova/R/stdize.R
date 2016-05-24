### stdize.R: center and scale a matrix.
### By Bjoern-Helge Mevik
### $Id: stdize.R 53 2007-04-20 12:05:00Z bhm $

stdize <- function(x, center = TRUE, scale = TRUE, avoid.zero.divisor = FALSE) {
    n <- nrow(x)
    ones <- matrix(1, nrow = n, ncol = 1)
    if (is.logical(center)) {
        if(center) {
            x <- x - ones %*% colMeans(x)
        }
    } else {
        x <- x - ones %*% center
    }
    if (is.logical(scale)) {
        if(scale) {
            varfun <-
                if(isTRUE(center)) function(v) sum(v^2) / (n - 1) else var
            scale <- sqrt(apply(x, 2, varfun))
            if (avoid.zero.divisor) scale[scale == 0] <- 1
            x <- x / (ones %*% scale)
        }
    } else {
        if (avoid.zero.divisor) scale[scale == 0] <- 1
        x <- x / (ones %*% scale)
    }
    x
}
