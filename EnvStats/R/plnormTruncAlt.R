plnormTruncAlt <-
function (q, mean = exp(1/2), cv = sqrt(exp(1) - 1), min = 0, 
    max = Inf) 
{
    names.q <- names(q)
    arg.mat <- cbind.no.warn(q = as.vector(q), mean = as.vector(mean), 
        cv = as.vector(cv), min = as.vector(min), max = as.vector(max))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        p <- rep(NA, nrow(arg.mat))
    else {
        p <- numeric(nrow(arg.mat))
        p[na.index] <- NA
        p.no.na <- p[!na.index]
        for (i in c("q", "mean", "cv", "min", "max")) assign(i, 
            arg.mat[!na.index, i])
        if (any(mean < .Machine$double.eps)) 
            stop("All non-missing values of 'mean' must be positive.")
        if (any(cv < .Machine$double.eps)) 
            stop("All non-missing values of 'cv' must be positive.")
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
            mean <- mean[index]
            cv <- cv[index]
            F.min <- plnormAlt(min[index], mean = mean, cv = cv)
            p.no.na[index] <- (plnormAlt(q[index], mean = mean, 
                cv = cv) - F.min)/(plnormAlt(max[index], mean = mean, 
                cv = cv) - F.min)
        }
        p[!na.index] <- p.no.na
    }
    if (!is.null(names.q)) 
        names(p) <- rep(names.q, length = length(p))
    p
}
