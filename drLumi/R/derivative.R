## Helper function to derivate of order 1, 2, 3 etc etc
## warning about approximations 
derivative <- function(f, x, ..., order = 1, delta = 0.1, sig = 6) {
    ## Numerically computes the specified order derivative of f at x
    vals <- matrix(NA, nrow = order + 1, ncol = order + 1)
    grid <- seq(x - delta/2, x + delta/2, length.out = order + 1)
    vals[1, ] <- sapply(grid, f, ...) - f(x, ...)
    for (i in 2:(order + 1)) {
        for (j in 1:(order - i + 2)) {
            stepsize <- grid[i + j - 1] - grid[i + j - 2]
            vals[i, j] <- (vals[i - 1, j + 1] - vals[i - 1, j])/stepsize
        }
    }
    return(signif(vals[order + 1, 1], sig))
}

