updateMean <- function (xbar, x, n, f, byrow = TRUE) 
{
    if (missing(n) && missing(f)) 
        stop("At least one of the arguments 'n' and 'f' must be specified")
    if (!is.matrix(x)) {
        if (missing(f)) 
        		f <- 1/(n+1)
        	return((1-f)*xbar + f*x)
    }
    dimx <- dim(x)
    k <- ifelse(byrow, dimx[1], dimx[2])
    p <- ifelse(byrow, dimx[2], dimx[1])
    if (length(xbar) != p) 
        stop("'x' and 'xbar' of incompatible dimensions.\nCheck these arguments and 'byrow'")
    if (missing(f)) 
        f <- k/(n + k)
    mean_x <- if (byrow) {
        .colMeans(x, k, p)
    }
    else .rowMeans(x, p, k)
    return((1 - f) * xbar + f * mean_x)
}