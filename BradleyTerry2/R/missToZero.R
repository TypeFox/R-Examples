missToZero <- function(x, miss, dim = 1) {
    if (dim == 1) x[miss, ] <- 0
    else x[, miss] <- 0
    x
}
