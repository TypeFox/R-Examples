pzmnorm <-
function (q, mean = 0, sd = 1, p.zero = 0.5) 
{
    names.q <- names(q)
    arg.mat <- cbind.no.warn(q = as.vector(q), mean = as.vector(mean), 
        sd = as.vector(sd), p.zero = as.vector(p.zero))
    for (i in c("q", "mean", "sd", "p.zero")) assign(i, arg.mat[, 
        i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        p <- rep(NA, length(q))
    else {
        if (any(sd[!na.index] < .Machine$double.eps)) 
            stop("All values of 'sd' must be positive.")
        if (any(p.zero[!na.index] <= 0) || any(p.zero[!na.index] >= 
            1)) 
            stop(paste("All values of 'p.zero' must be", "greater than 0 and less than 1."))
        p <- (1 - p.zero) * pnorm(q, mean, sd)
        index <- !na.index & q >= 0
        p[index] <- p[index] + p.zero[index]
    }
    if (!is.null(names.q)) 
        names(p) <- rep(names.q, length = length(p))
    else names(p) <- NULL
    p
}
