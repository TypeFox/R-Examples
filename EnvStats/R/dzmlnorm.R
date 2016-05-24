dzmlnorm <-
function (x, meanlog = 0, sdlog = 1, p.zero = 0.5) 
{
    names.x <- names(x)
    arg.mat <- cbind.no.warn(x = as.vector(x), meanlog = as.vector(meanlog), 
        sdlog = as.vector(sdlog), p.zero = as.vector(p.zero))
    for (i in c("x", "meanlog", "sdlog", "p.zero")) assign(i, 
        arg.mat[, i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        y <- rep(NA, length(x))
    else {
        if (any(sdlog[!na.index] < .Machine$double.eps)) 
            stop("All values of 'sdlog' must be positive.")
        if (any(p.zero[!na.index] <= 0) || any(p.zero[!na.index] >= 
            1)) 
            stop(paste("All values of 'p.zero' must be", "greater than 0 and less than 1."))
        y <- (1 - p.zero) * dlnorm(x, meanlog, sdlog)
        index <- !na.index & x == 0
        y[index] <- p.zero[index]
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(y))
    else names(y) <- NULL
    y
}
