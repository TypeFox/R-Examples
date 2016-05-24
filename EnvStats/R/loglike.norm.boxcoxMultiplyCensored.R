loglike.norm.boxcoxMultiplyCensored <-
function (x, censored, censoring.side, lambda, mean, sd, eps = .Machine$double.eps) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    if (!(is.vector(censored, mode = "numeric") || is.vector(censored, 
        mode = "logical"))) 
        stop("'censored' must be a logical or numeric vector")
    if (length(censored) != length(x)) 
        stop("'censored' must be the same length as 'x'")
    if (any(is.na(censored))) 
        stop("'censored' cannot contain missing values")
    if (is.numeric(censored)) {
        if (!all(censored == 0 | censored == 1)) 
            stop(paste("When 'censored' is a numeric vector, all values of", 
                "'censored' must be 0 (not censored) or 1 (censored)."))
        censored <- as.logical(censored)
    }
    n.cen <- sum(censored)
    if (n.cen == 0) 
        stop("No censored values indicated by 'censored'.")
    censoring.side <- match.arg(censoring.side, c("left", "right"))
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
    censoring.name <- deparse(substitute(censored))
    if (any(is.na(x))) {
        statistic <- NA
    }
    else {
        if (any(x <= 0)) 
            stop("All values of 'x' must be positive")
        x.no.cen <- x[!censored]
        if (length(unique(x.no.cen)) < 2) 
            stop(paste("'x' must contain at least 2 non-missing,", 
                "uncensored, distinct values."))
        x.cen <- x[censored]
        cen.levels <- sort(unique(x.cen))
        N <- length(x)
        y <- boxcoxTransform(x = x, lambda = lambda, eps = eps)
        statistic <- loglike.norm.multiply.censored(x = y, censored = censored, 
            censoring.side = censoring.side, mean = mean, sd = sd)$statistic + 
            (lambda - 1) * sum(log(x.no.cen))
    }
    statistic
}
