##
##  n c h o o s e k . R  Binomial Coefficients
##


nchoosek <- function(n, k) {
    stopifnot(is.numeric(n), length(n) == 1, n >= 0, floor(n) == ceiling(n),
              is.numeric(k), length(k) == 1, k >= 0, floor(k) == ceiling(k))
    if (k > n)
        stop("Argument 'k' must be an integer between 0 and 'n'.")

    choose(n, k)
}
