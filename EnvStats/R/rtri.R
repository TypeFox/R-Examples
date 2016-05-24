rtri <-
function (n, min = 0, max = 1, mode = 1/2) 
{
    ln <- length(n)
    if (ln < 1) 
        stop("'n' must be non-empty.")
    if (ln > 1) 
        n <- ln
    else {
        if (is.na(n) || n <= 0 || n != trunc(n)) 
            stop("'n' must be a positive integer or vector.")
    }
    arg.mat <- cbind.no.warn(dum = rep(1, n), min = as.vector(min), 
        max = as.vector(max), mode = as.vector(mode))[, -1, drop = FALSE]
    if (n < nrow(arg.mat)) 
        arg.mat <- arg.mat[1:n, , drop = FALSE]
    for (i in c("min", "max", "mode")) assign(i, arg.mat[, i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        return(rep(NA, n))
    else {
        if (any(is.infinite(min)) || any(is.infinite(max))) 
            stop("All non-missing values of 'min' and 'max' must be finite.")
        if (any(mode <= min) || any(max <= mode)) 
            stop(paste("All values of 'mode' must be larger than", 
                "the corresponding values of 'min', and all", 
                "values of 'max' must be larger than the", "corresponding values of 'mode'."))
        return(qtri(p = runif(n), min = min, max = max, mode = mode))
    }
}
