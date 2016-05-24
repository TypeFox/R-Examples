dlnormTrunc <-
function (x, meanlog = 0, sdlog = 1, min = 0, max = Inf) 
{
    names.x <- names(x)
    arg.mat <- cbind.no.warn(x = as.vector(x), meanlog = as.vector(meanlog), 
        sdlog = as.vector(sdlog), min = as.vector(min), max = as.vector(max))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        y <- rep(NA, nrow(arg.mat))
    else {
        y <- numeric(nrow(arg.mat))
        y[na.index] <- NA
        y.no.na <- y[!na.index]
        for (i in c("x", "meanlog", "sdlog", "min", "max")) assign(i, 
            arg.mat[!na.index, i])
        if (any(sdlog < .Machine$double.eps)) 
            stop("All non-missing values of 'sdlog' must be positive.")
        if (any(min >= max)) 
            stop(paste("All non-missing values of 'min' must be", 
                "less than the corresponding elements of 'max'."))
        if (any(min < 0)) 
            stop(paste("All non-missing values of 'min' must be", 
                "greater than or equal to 0."))
        x.out <- x < min | x > max
        y.no.na[x.out] <- 0
        if (any(index <- !x.out)) {
            meanlog <- meanlog[index]
            sdlog <- sdlog[index]
            y.no.na[index] <- dlnorm(x[index], meanlog = meanlog, 
                sdlog = sdlog)/(plnorm(max[index], meanlog = meanlog, 
                sdlog = sdlog) - plnorm(min[index], meanlog = meanlog, 
                sdlog = sdlog))
        }
        y[!na.index] <- y.no.na
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(y))
    y
}
