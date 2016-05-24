##
##  d e t r e n d . R  Remove Linear Trends
##


detrend <- function(x, tt = 'linear', bp = c()) {
    if (!is.numeric(x) && !is.complex(x))
        stop("'x' must be a numeric or complex vector or matrix.")
    trendType <- pmatch(tt, c('constant', 'linear'), nomatch = 0)

    if (is.vector(x))
        x <- as.matrix(x)
    n <- nrow(x)
    if (length(bp) > 0 && !all(bp %in% 1:n))
        stop("Breakpoints 'bp' must elements of 1:length(x).")

    if (trendType == 1) {  # 'constant'
        if (!is.null(bp))
            warning("Breakpoints not used for 'constant' trend type.")
        y <- x - matrix(1, n, 1) %*% apply(x, 2, mean)

    } else if (trendType == 2) {  # 'linear'
        bp <- sort(unique(c(0, c(bp), n-1)))
        lb <- length(bp) - 1

        a <- cbind(matrix(0, n, lb), matrix(1, n, 1))
        for (kb in 1:lb) {
            m <- n - bp[kb]
            a[(1:m) + bp[kb], kb] <- as.matrix(1:m)/m
        }
        y <- x - a %*% qr.solve(a, x)

    } else {
        stop("Trend type 'tt' must be 'constant' or 'linear'.")
    }
    
    return(y)
}
