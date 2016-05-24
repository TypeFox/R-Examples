dzmnorm <-
function (x, mean = 0, sd = 1, p.zero = 0.5) 
{
    names.x <- names(x)
    arg.mat <- cbind.no.warn(x = as.vector(x), mean = as.vector(mean), 
        sd = as.vector(sd), p.zero = as.vector(p.zero))
    for (i in c("x", "mean", "sd", "p.zero")) assign(i, arg.mat[, 
        i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        y <- rep(NA, length(x))
    else {
        if (any(sd[!na.index] < .Machine$double.eps)) 
            stop("All values of 'sd' must be positive.")
        if (any(p.zero[!na.index] <= 0) || any(p.zero[!na.index] >= 
            1)) 
            stop(paste("All values of 'p.zero' must be", "greater than 0 and less than 1."))
        y <- (1 - p.zero) * dnorm(x, mean, sd)
        y[!na.index & x == 0] <- p.zero[!na.index & x == 0]
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(y))
    else names(y) <- NULL
    y
}
