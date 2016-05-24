rlnormAlt <-
function (n, mean = exp(1/2), cv = sqrt(exp(1) - 1)) 
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
    arg.mat <- cbind.no.warn(dum = rep(1, n), mean = as.vector(mean), 
        cv = as.vector(cv))[, -1, drop = FALSE]
    if (n < nrow(arg.mat)) 
        arg.mat <- arg.mat[1:n, , drop = FALSE]
    for (i in c("mean", "cv")) assign(i, arg.mat[, i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        return(rep(NA, n))
    else {
        if (any(c(mean[!na.index], cv[!na.index]) < .Machine$double.eps)) 
            stop("All values of 'mean' and 'cv' must be positive.")
        sdlog <- sqrt(log(1 + cv^2))
        meanlog <- log(mean) - (sdlog^2)/2
        return(rlnorm(n, meanlog, sdlog))
    }
}
