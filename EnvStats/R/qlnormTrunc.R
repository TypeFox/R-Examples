qlnormTrunc <-
function (p, meanlog = 0, sdlog = 1, min = 0, max = Inf) 
{
    names.p <- names(p)
    arg.mat <- cbind.no.warn(p = as.vector(p), meanlog = as.vector(meanlog), 
        sdlog = as.vector(sdlog), min = as.vector(min), max = as.vector(max))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        q <- rep(NA, nrow(arg.mat))
    else {
        q <- numeric(nrow(arg.mat))
        q[na.index] <- NA
        q.no.na <- q[!na.index]
        for (i in c("p", "meanlog", "sdlog", "min", "max")) assign(i, 
            arg.mat[!na.index, i])
        if (any(p < 0) || any(p > 1)) 
            stop("All non-missing values of 'p' must be between 0 and 1")
        if (any(sdlog < .Machine$double.eps)) 
            stop("All non-missing values of 'sdlog' must be positive.")
        if (any(min >= max)) 
            stop(paste("All non-missing values of 'min' must be", 
                "less than the corresponding elements of 'max'."))
        if (any(min < 0)) 
            stop(paste("All non-missing values of 'min' must be", 
                "greater than or equal to 0."))
        p.low <- p == 0
        q.no.na[p.low] <- min[p.low]
        p.high <- p == 1
        q.no.na[p.high] <- max[p.high]
        if (any(index <- !(p.low | p.high))) {
            meanlog <- meanlog[index]
            sdlog <- sdlog[index]
            F.min <- plnorm(min[index], meanlog = meanlog, sdlog = sdlog)
            q.no.na[index] <- qlnorm(p[index] * (plnorm(max[index], 
                meanlog = meanlog, sdlog = sdlog) - F.min) + 
                F.min, meanlog = meanlog, sdlog = sdlog)
        }
        q[!na.index] <- q.no.na
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    q
}
