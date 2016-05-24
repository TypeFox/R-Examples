boxcoxTransform <-
function (x, lambda, eps = .Machine$double.eps) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    if (!is.vector(lambda, mode = "numeric") || length(lambda) != 
        1 || !is.finite(lambda)) 
        stop("'lambda' must be a non-missing, finite numeric scalar")
    if (!is.vector(eps, mode = "numeric") || length(eps) != 1 || 
        !is.finite(eps) || eps <= 0) 
        stop("'eps' must be a non-missing, finite, positive numeric scalar")
    if (all(is.na(x))) 
        new.x <- rep(NA, length(x))
    else {
        if (any(x[!is.na(x)] <= 0)) 
            stop("All non-missing values of 'x' must be positive")
        if (abs(lambda) < eps) 
            y <- log(x)
        else y <- (x^lambda - 1)/lambda
    }
    y
}
