rpareto <-
function (n, location, shape = 1) 
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
    arg.mat <- cbind.no.warn(dum = rep(1, n), location = as.vector(location), 
        shape = as.vector(shape))[, -1, drop = FALSE]
    if (n < nrow(arg.mat)) 
        arg.mat <- arg.mat[1:n, , drop = FALSE]
    for (i in c("location", "shape")) assign(i, arg.mat[, i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        return(rep(NA, n))
    else {
        if (any(c(location[!na.index], shape[!na.index]) < .Machine$double.eps)) 
            stop("All values of 'location' and 'shape' must be positive.")
        return(qpareto(p = runif(n), location = location, shape = shape))
    }
}
