rgevd <-
function (n, location = 0, scale = 1, shape = 0) 
{
    ln <- length(n)
    if (ln < 1) 
        stop("'n' must be non-empty")
    if (ln > 1) 
        n <- ln
    else {
        if (is.na(n) || n <= 0 || n != trunc(n)) 
            stop("'n' must be a positive integer or a vector.")
    }
    arg.mat <- cbind.no.warn(dum = rep(1, n), location = as.vector(location), 
        scale = as.vector(scale), shape = as.vector(shape))[, 
        -1, drop = FALSE]
    if (n < nrow(arg.mat)) 
        arg.mat <- arg.mat[1:n, , drop = FALSE]
    for (i in c("location", "scale", "shape")) assign(i, arg.mat[, 
        i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        return(rep(NA, n))
    else {
        if (any(scale[!na.index] < .Machine$double.eps)) 
            stop("All values of 'scale' must be positive.")
        return(qgevd(p = runif(n), location = location, scale = scale, 
            shape = shape))
    }
}
