dzmlnormAlt <-
function (x, mean = exp(1/2), cv = sqrt(exp(1) - 1), p.zero = 0.5) 
{
    names.x <- names(x)
    arg.mat <- cbind.no.warn(x = as.vector(x), mean = as.vector(mean), 
        cv = as.vector(cv), p.zero = as.vector(p.zero))
    for (i in c("x", "mean", "cv", "p.zero")) assign(i, arg.mat[, 
        i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        y <- rep(NA, length(x))
    else {
        if (any(c(mean[!na.index], cv[!na.index]) < .Machine$double.eps)) 
            stop("All values of 'mean' and 'cv' must be positive.")
        if (any(p.zero[!na.index] <= 0) || any(p.zero[!na.index] >= 
            1)) 
            stop(paste("All non-missing values of 'p.zero' must be", 
                "greater than 0 and less than 1."))
        y <- (1 - p.zero) * dlnormAlt(x, mean, cv)
        index <- !na.index & x == 0
        y[index] <- p.zero[index]
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(y))
    else names(y) <- NULL
    y
}
