loglike.norm.boxcox <-
function (x, lambda, mean, sd, eps = .Machine$double.eps) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    if (!is.vector(lambda, mode = "numeric") || length(lambda) != 
        1 || !is.finite(lambda)) 
        stop("'lambda' must be a non-missing, finite numeric scalar")
    if (!is.vector(mean, mode = "numeric") || length(mean) != 
        1 || !is.finite(mean)) 
        stop("'mean' must be a non-missing, finite numeric scalar")
    if (!is.vector(sd, mode = "numeric") || length(sd) != 1 || 
        !is.finite(sd) || sd < 0) 
        stop("'sd' must be a non-missing, finite, positive numeric scalar")
    data.name <- deparse(substitute(x))
    n <- length(x)
    if (any(is.na(x))) {
        statistic <- NA
    }
    else {
        if (any(x <= 0)) 
            stop("All values of 'x' must be positive")
        y <- boxcoxTransform(x = x, lambda = lambda, eps = eps)
        statistic <- loglike.norm(x = y, mean = mean, sd = sd)$statistic + 
            (lambda - 1) * sum(log(x))
    }
    statistic
}
