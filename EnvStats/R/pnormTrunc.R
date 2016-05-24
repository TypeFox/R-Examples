pnormTrunc <-
function (q, mean = 0, sd = 1, min = -Inf, max = Inf) 
{
    names.q <- names(q)
    arg.mat <- cbind.no.warn(q = as.vector(q), mean = as.vector(mean), 
        sd = as.vector(sd), min = as.vector(min), max = as.vector(max))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        p <- rep(NA, nrow(arg.mat))
    else {
        p <- numeric(nrow(arg.mat))
        p[na.index] <- NA
        p.no.na <- p[!na.index]
        for (i in c("q", "mean", "sd", "min", "max")) assign(i, 
            arg.mat[!na.index, i])
        if (any(sd < .Machine$double.eps)) 
            stop("All non-missing values of 'sd' must be positive.")
        if (any(min >= max)) 
            stop(paste("All non-missing values of 'min' must be", 
                "less than the corresponding elements of 'max'."))
        q.low <- q < min
        p.no.na[q.low] <- 0
        q.high <- q > max
        p.no.na[q.high] <- 1
        if (any(index <- !(q.low | q.high))) {
            mean <- mean[index]
            sd <- sd[index]
            F.min <- pnorm(min[index], mean = mean, sd = sd)
            p.no.na[index] <- (pnorm(q[index], mean = mean, sd = sd) - 
                F.min)/(pnorm(max[index], mean = mean, sd = sd) - 
                F.min)
        }
        p[!na.index] <- p.no.na
    }
    if (!is.null(names.q)) 
        names(p) <- rep(names.q, length = length(p))
    p
}
