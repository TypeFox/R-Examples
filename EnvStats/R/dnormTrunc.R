dnormTrunc <-
function (x, mean = 0, sd = 1, min = -Inf, max = Inf) 
{
    names.x <- names(x)
    arg.mat <- cbind.no.warn(x = as.vector(x), mean = as.vector(mean), 
        sd = as.vector(sd), min = as.vector(min), max = as.vector(max))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        y <- rep(NA, nrow(arg.mat))
    else {
        y <- numeric(nrow(arg.mat))
        y[na.index] <- NA
        y.no.na <- y[!na.index]
        for (i in c("x", "mean", "sd", "min", "max")) assign(i, 
            arg.mat[!na.index, i])
        if (any(sd < .Machine$double.eps)) 
            stop("All non-missing values of 'sd' must be positive.")
        if (any(min >= max)) 
            stop(paste("All non-missing values of 'min' must be", 
                "less than the corresponding elements of 'max'."))
        x.out <- x < min | x > max
        y.no.na[x.out] <- 0
        if (any(index <- !x.out)) {
            mean <- mean[index]
            sd <- sd[index]
            y.no.na[index] <- dnorm(x[index], mean = mean, sd = sd)/(pnorm(max[index], 
                mean = mean, sd = sd) - pnorm(min[index], mean = mean, 
                sd = sd))
        }
        y[!na.index] <- y.no.na
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(y))
    y
}
