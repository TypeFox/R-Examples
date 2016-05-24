dlnormAlt <-
function (x, mean = exp(1/2), cv = sqrt(exp(1) - 1), log = FALSE) 
{
    names.x <- names(x)
    arg.mat <- cbind.no.warn(x = as.vector(x), mean = as.vector(mean), 
        cv = as.vector(cv))
    for (i in c("x", "mean", "cv")) assign(i, arg.mat[, i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        y <- rep(NA, length(x))
    else {
        if (any(c(mean[!na.index], cv[!na.index]) < .Machine$double.eps)) 
            stop("All values of 'mean' and 'cv' must be positive.")
        sdlog <- sqrt(log(1 + cv^2))
        meanlog <- log(mean) - (sdlog^2)/2
        y <- dlnorm(x = x, meanlog = meanlog, sdlog = sdlog, 
            log = log)
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(y))
    else names(y) <- NULL
    y
}
