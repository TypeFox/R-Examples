qlnormTruncAlt <-
function (p, mean = exp(1/2), cv = sqrt(exp(1) - 1), min = 0, 
    max = Inf) 
{
    names.p <- names(p)
    arg.mat <- cbind.no.warn(p = as.vector(p), mean = as.vector(mean), 
        cv = as.vector(cv), min = as.vector(min), max = as.vector(max))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        q <- rep(NA, nrow(arg.mat))
    else {
        q <- numeric(nrow(arg.mat))
        q[na.index] <- NA
        q.no.na <- q[!na.index]
        for (i in c("p", "mean", "cv", "min", "max")) assign(i, 
            arg.mat[!na.index, i])
        if (any(p < 0) || any(p > 1)) 
            stop("All non-missing values of 'p' must be between 0 and 1")
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
        p.low <- p == 0
        q.no.na[p.low] <- min[p.low]
        p.high <- p == 1
        q.no.na[p.high] <- max[p.high]
        if (any(index <- !(p.low | p.high))) {
            mean <- mean[index]
            cv <- cv[index]
            F.min <- plnormAlt(min[index], mean = mean, cv = cv)
            q.no.na[index] <- qlnormAlt(p[index] * (plnormAlt(max[index], 
                mean = mean, cv = cv) - F.min) + F.min, mean = mean, 
                cv = cv)
        }
        q[!na.index] <- q.no.na
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    q
}
