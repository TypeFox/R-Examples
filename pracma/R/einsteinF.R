##
##  e i n s t e i n F . R  Einstein Functions
##


einsteinF <- function(d, x) {
    stopifnot(is.numeric(x) || is.complex(x))
    fi <- which(x == 0)

    if (d == 1) {
        y <- x^2 * exp(x) / (exp(x) - 1)^2
        y[fi] <- 1
    } else if (d == 2) {
        y <- x / (exp(x) - 1)
        y[fi] <- 1
    } else if (d == 3) {
        y <- log(1 - exp(-x))
        y[fi] <- -Inf
    } else if (d == 4) {
        y <- x / (exp(x) - 1) - log(1 - exp(-x))
        y[fi] <- Inf
    } else {
        stop("Argument 'd' must be one of the integers 1, 2, 3, 4.")
    }

    return(y)
}
