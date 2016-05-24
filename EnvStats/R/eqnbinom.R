eqnbinom <-
function (x, size = NULL, p = 0.5, method = "mle/mme", digits = 0) 
{
    if (!is.vector(p, mode = "numeric") || is.factor(p)) 
        stop("'p' must be a numeric vector.")
    if (any(!is.finite(p))) 
        stop("NA/NaN/Inf values not allowed in 'p'.")
    if (any(p < 0) || any(p > 1)) 
        stop("All values of 'p' must be between 0 and 1.")
    method <- match.arg(method, c("mle/mme", "mvue"))
    if (x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored") {
        if (x$distribution != "Negative Binomial") 
            stop(paste("'eqnbinom' estimates quantiles", "for a negative binomial distribution.  You have supplied an object", 
                "that assumes a different distribution."))
        class.x <- oldClass(x)
        if (!is.null(x$interval)) {
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        size <- x$parameters["size"]
        prob <- x$parameters["prob"]
        n <- x$sample.size
        ret.list <- x
    }
    else {
        if (!is.vector(x, mode = "numeric") || !is.vector(size, 
            mode = "numeric")) 
            stop(paste("'x' must be either a list that inherits from", 
                "the class 'estimate', or else a numeric vector,", 
                "and 'size' must be a numeric vector"))
        data.name <- deparse(substitute(x))
        if ((bad.obs <- sum(!(all.ok <- is.finite(x) & is.finite(size)))) > 
            0) {
            is.not.finite.warning(x)
            is.not.finite.warning(size)
            x <- x[all.ok]
            size <- size[all.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' and/or 'size' removed."))
        }
        n <- length(x)
        if (n < 1) 
            stop("'x' and 'size' must contain at least one non-missing pair of values.")
        if (!all(x == trunc(x)) || any(x < 0) || !all(size == 
            trunc(size)) || any(size < 1)) 
            stop(paste("All values of 'x' must be non-negative integers,", 
                "and all values of 'size' must be positive integers."))
        ret.list <- enbinom(x, size = size, method = method)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        size <- ret.list$parameters["size"]
        prob <- ret.list$parameters["prob"]
    }
    q <- qnbinom(p, size = size, prob = prob)
    if (length(p) == 1 && p == 0.5) 
        names(q) <- "Median"
    else {
        pct <- round(100 * p, digits)
        names(q) <- paste(pct, number.suffix(pct), " %ile", sep = "")
    }
    ret.list <- c(ret.list, list(quantiles = q))
    ret.list$quantile.method <- paste("Quantile(s) Based on\n", 
        space(33), ret.list$method, " Estimators", sep = "")
    if (x.is.est.obj) 
        oldClass(ret.list) <- class.x
    else oldClass(ret.list) <- "estimate"
    ret.list
}
