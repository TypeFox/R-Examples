###
### NTHROOT.R - Compute the real n-th root
###


nthroot <- function(x, n) {
    if (! is.numeric(x))
        stop("Argument 'x' must be numeric.")
    if (missing(n) || n <= 0 || ceiling(n) != floor(n))
        stop("Argument 'n' must be a positive integer.")
    if (any(x[!is.na(x)] < 0) && n %% 2 == 0)
        stop("If argument 'x' is negative, 'n' must be an odd integer.")

    sx <- sign(x)
    return(sx * (sx * x)^(1/n))
}
