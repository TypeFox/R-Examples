hatMat <- function(x, trace = FALSE,
                   pred.sm = function(x,y,...)
                   predict(smooth.spline(x,y, ...), x = x)$y,
                   ...)
{
    ## Purpose: Return Hat matrix of a smoother -- very general (but slow)
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  7 Mar 2001, 11:12
    stopifnot(is.logical(trace), length(trace) == 1)
    n <- NROW(x)
    if(is.unsorted(x) && !missing(pred.sm))
	warning("'x' is not sorted increasingly:\n ",
		"  this may be inefficient and lead to wrong results")
    y <- pred.sm(x, numeric(n), ...)
    if(!is.numeric(y) || length(y) !=n)
        stop("`pred.sm' does not return a numeric length n vector")
    H <- if(trace) 0 else matrix(as.numeric(NA), n,n)
    for (i in 1:n) {
        y <- numeric(n) ; y[i] <- 1
        y <- pred.sm(x, y, ...)
        if(trace) H <- H + y[i] else H[,i] <- y
    }
    H
}

