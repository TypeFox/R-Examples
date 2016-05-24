plnormTrunc <-
function (q, meanlog = 0, sdlog = 1, min = 0, max = Inf) 
{
    names.q <- names(q)
    arg.mat <- cbind.no.warn(q = as.vector(q), meanlog = as.vector(meanlog), 
        sdlog = as.vector(sdlog), min = as.vector(min), max = as.vector(max))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        p <- rep(NA, nrow(arg.mat))
    else {
        p <- numeric(nrow(arg.mat))
        p[na.index] <- NA
        p.no.na <- p[!na.index]
        for (i in c("q", "meanlog", "sdlog", "min", "max")) assign(i, 
            arg.mat[!na.index, i])
        if (any(sdlog < .Machine$double.eps)) 
            stop("All non-missing values of 'sdlog' must be positive.")
        if (any(min >= max)) 
            stop(paste("All non-missing values of 'min' must be", 
                "less than the corresponding elements of 'max'."))
        if (any(min < 0)) 
            stop(paste("All non-missing values of 'min' must be", 
                "greater than or equal to 0."))
        q.low <- q < min
        p.no.na[q.low] <- 0
        q.high <- q > max
        p.no.na[q.high] <- 1
        if (any(index <- !(q.low | q.high))) {
            meanlog <- meanlog[index]
            sdlog <- sdlog[index]
            F.min <- plnorm(min[index], meanlog = meanlog, sdlog = sdlog)
            p.no.na[index] <- (plnorm(q[index], meanlog = meanlog, 
                sdlog = sdlog) - F.min)/(plnorm(max[index], meanlog = meanlog, 
                sdlog = sdlog) - F.min)
        }
        p[!na.index] <- p.no.na
    }
    if (!is.null(names.q)) 
        names(p) <- rep(names.q, length = length(p))
    p
}
