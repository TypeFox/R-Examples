eqzmlnormAlt <-
function (x, p = 0.5, method = "mvue", digits = 0) 
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
        if (x$distribution != "Zero-Modified Lognormal (Delta)" || 
            names(x$parameters) != c("mean", "cv", "p.zero", 
                "mean.zmlnorm", "cv.zmlnorm")) 
            stop(paste("'eqzmlnormAlt' estimates quantiles", 
                "for a zero-modified lognormal distribution", 
                "(alternative parameterization).", "You have supplied an object", 
                "that assumes a different distribution."))
        class.x <- oldClass(x)
        if (!is.null(x$interval)) {
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        mean <- x$parameters["mean"]
        cv <- x$parameters["cv"]
        p.zero <- x$parameters["p.zero"]
        n <- x$sample.size
        ret.list <- x
    }
    else {
        if (!is.vector(x, mode = "numeric")) 
            stop(paste("'x' must be either a list that inherits from", 
                "the class 'estimate', or else a numeric vector"))
        data.name <- deparse(substitute(x))
        if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
            is.not.finite.warning(x)
            x <- x[x.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
        }
        n <- length(x)
        if (n < 1 || any(x < 0)) 
            stop(paste("'x' must contain at least 1 non-missing value,", 
                "and all non-missing values of 'x' must be non-negative. ", 
                "This is not true for 'x' =", data.name))
        ret.list <- ezmlnormAlt(x, method = method)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        mean <- ret.list$parameters["mean"]
        cv <- ret.list$parameters["cv"]
        p.zero <- ret.list$parameters["p.zero"]
    }
    q <- qzmlnormAlt(p, mean = mean, cv = cv, p.zero = p.zero)
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
