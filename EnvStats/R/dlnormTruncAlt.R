dlnormTruncAlt <-
function (x, mean = exp(1/2), cv = sqrt(exp(1) - 1), min = 0, 
    max = Inf) 
{
    names.x <- names(x)
    arg.mat <- cbind.no.warn(x = as.vector(x), mean = as.vector(mean), 
        cv = as.vector(cv), min = as.vector(min), max = as.vector(max))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        y <- rep(NA, nrow(arg.mat))
    else {
        y <- numeric(nrow(arg.mat))
        y[na.index] <- NA
        y.no.na <- y[!na.index]
        for (i in c("x", "mean", "cv", "min", "max")) assign(i, 
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
        x.out <- x < min | x > max
        y.no.na[x.out] <- 0
        if (any(index <- !x.out)) {
            mean <- mean[index]
            cv <- cv[index]
            y.no.na[index] <- dlnormAlt(x[index], mean = mean, 
                cv = cv)/(plnormAlt(max[index], mean = mean, 
                cv = cv) - plnormAlt(min[index], mean = mean, 
                cv = cv))
        }
        y[!na.index] <- y.no.na
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(y))
    y
}
