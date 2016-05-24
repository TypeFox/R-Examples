qnormTrunc <-
function (p, mean = 0, sd = 1, min = -Inf, max = Inf) 
{
    names.p <- names(p)
    arg.mat <- cbind.no.warn(p = as.vector(p), mean = as.vector(mean), 
        sd = as.vector(sd), min = as.vector(min), max = as.vector(max))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        q <- rep(NA, nrow(arg.mat))
    else {
        q <- numeric(nrow(arg.mat))
        q[na.index] <- NA
        q.no.na <- q[!na.index]
        for (i in c("p", "mean", "sd", "min", "max")) assign(i, 
            arg.mat[!na.index, i])
        if (any(p < 0) || any(p > 1)) 
            stop("All non-missing values of 'p' must be between 0 and 1")
        if (any(sd < .Machine$double.eps)) 
            stop("All non-missing values of 'sd' must be positive.")
        if (any(min >= max)) 
            stop(paste("All non-missing values of 'min' must be", 
                "less than the corresponding elements of 'max'."))
        p.low <- p == 0
        q.no.na[p.low] <- min[p.low]
        p.high <- p == 1
        q.no.na[p.high] <- max[p.high]
        if (any(index <- !(p.low | p.high))) {
            mean <- mean[index]
            sd <- sd[index]
            F.min <- pnorm(min[index], mean = mean, sd = sd)
            q.no.na[index] <- qnorm(p[index] * (pnorm(max[index], 
                mean = mean, sd = sd) - F.min) + F.min, mean = mean, 
                sd = sd)
        }
        q[!na.index] <- q.no.na
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    q
}
