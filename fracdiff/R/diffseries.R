#### by Valderio Reisen  -- Dec.2005

## MM:  This is 'not optimal' -- and I may have better in ../filters.R ? <<< FIXME >>>
diffseries <- function(x, d)
{
    x <- as.data.frame(x)
    names(x) <- "series"
    x <- x$series
    if (NCOL(x) > 1)
        stop("only implemented for univariate time series")
    if (any(is.na(x)))
        stop("NAs in x")
    n <- length(x)
    stopifnot(n >= 2)
    x <- x - mean(x)
    PI <- numeric(n)
    PI[1] <- -d
    for (k in 2:n) {
        PI[k] <- PI[k-1]*(k - 1 - d)/k
    }
    ydiff <- x
    for (i in 2:n) {
        ydiff[i] <- x[i] + sum(PI[1:(i-1)]*x[(i-1):1])
    }
    ## return numeric!
    ydiff
}

