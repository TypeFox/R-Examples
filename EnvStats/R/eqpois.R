eqpois <-
function (x, p = 0.5, method = "mle/mme/mvue", ci = FALSE, ci.method = "exact", 
    ci.type = "two-sided", conf.level = 0.95, digits = 0) 
{
    if (!is.vector(p, mode = "numeric") || is.factor(p)) 
        stop("'p' must be a numeric vector.")
    if (any(!is.finite(p))) 
        stop("NA/NaN/Inf values not allowed in 'p'.")
    if (any(p < 0) || any(p > 1)) 
        stop("All values of 'p' must be between 0 and 1.")
    method <- match.arg(method)
    if (x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored") {
        if (ci) 
            stop(paste("When ci=T you must supply", "the original observations; you cannot supply the", 
                "result of calling 'epois'"))
        if (x$distribution != "Poisson") 
            stop(paste("'eqpois' estimates quantiles", "for a Poisson distribution.  You have supplied an object", 
                "that assumes a different distribution."))
        class.x <- oldClass(x)
        if (!is.null(x$interval)) {
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        ret.list <- x
    }
    else {
        if (!is.vector(x, mode = "numeric") || is.factor(x)) 
            stop(paste("'x' must be either a list that inherits from", 
                "the class 'estimate', or else a numeric vector"))
        data.name <- deparse(substitute(x))
        if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
            is.not.finite.warning(x)
            x <- x[x.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
        }
        n <- length(x)
        if (n < 2 || length(unique(x)) < 2) 
            stop(paste("'x' must contain at least 2 non-missing distinct values. ", 
                "This is not true for 'x' =", data.name))
        if (any(x < 0) || any(x != trunc(x))) 
            stop("All non-missing values of 'x' must be non-negative integers")
        if (ci) {
            ci.type <- match.arg(ci.type, c("two-sided", "lower", 
                "upper"))
            ret.list <- epois(x, method = method, ci = TRUE, 
                ci.type = ci.type, conf.level = conf.level)
            conf.limits <- ret.list$interval$limits
            ret.list <- ret.list[-length(ret.list)]
        }
        else ret.list <- epois(x, method = method)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
    }
    lambda.hat <- ret.list$parameters
    q <- qpois(p, lambda = lambda.hat)
    if (length(p) == 1 && p == 0.5) 
        names(q) <- "Median"
    else {
        pct <- round(100 * p, digits)
        names(q) <- paste(pct, number.suffix(pct), " %ile", sep = "")
    }
    quantile.method <- "mle"
    ret.list <- c(ret.list, list(quantiles = q, quantile.method = quantile.method))
    if (ci) {
        if (length(p) > 1 || p <= 0 || p >= 1) 
            stop(paste("When 'ci' = TRUE, 'p' must be a scalar", 
                "larger than 0 and less than 1."))
        ci.method <- match.arg(ci.method)
        ci.type <- match.arg(ci.type)
        if (conf.level <= 0 || conf.level >= 1) 
            stop("The value of 'conf.level' must be between 0 and 1.")
        ci.obj <- ci.qpois(p, lambda.hat = lambda.hat, conf.limits = conf.limits, 
            n, ci.type, alpha = 1 - conf.level, digits)
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    if (x.is.est.obj) 
        oldClass(ret.list) <- class.x
    else oldClass(ret.list) <- "estimate"
    ret.list
}
