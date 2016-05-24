mknapsack <- function(p, w, k, bck = -1) {
    stopifnot(is.numeric(p), is.numeric(w), is.numeric(k))
    if (any(w <= 0))
        stop("'weights' must be a vector of positive numbers.")
    if (any(p <= 0))
        stop("'profits' must be a vector of positive numbers.")

    if (any(floor(p) != ceiling(p)) ||
        any(floor(w) != ceiling(w)) ||
        any(floor(k) != ceiling(k)) ||
        any(p >= 2^31) || any(w >= 2^31) || any(k >= 2^31))
        stop("All inputs must be positive integers < 2^31 !")

    n <- length(p); m <- length(k)
    if (length(w) != n)
        stop("Profit 'p' and weight 'w' must be vectors of equal length.")
    # bck <- -1, i.e., request exact solution,
    # else maximal number of allowed backtrackings

    xstar <- vector("integer", n)   # no. of knapsack item 'i' is assigned to
    vstar <- 0                      # value of optimal solution
    
    # dummy variables
    num <- 5*m + 14*n + 4*m*n + 3
    wk  <- numeric(n)
    iwk <- vector("integer", num)

    S <- .Fortran("mkp", as.integer(n), as.integer(m),
                    as.integer(p), as.integer(w), as.integer(k),
                    bs = as.integer(bck),
                    xs = as.integer(xstar), vs = as.integer(vstar),
                    as.numeric(wk), as.integer(iwk), as.integer(num),
                    PACKAGE = 'adagio')

    if (S$vs < 0)
        warning("Error condition raised: check input data ...!")

    return(list(ksack = S$xs, value = S$vs, btracks = S$bs))
}
