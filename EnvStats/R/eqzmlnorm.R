eqzmlnorm <-
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
            names(x$parameters) != c("meanlog", "sdlog", "p.zero", 
                "mean.zmlnorm", "sd.zmlnorm")) 
            stop(paste("'eqzmlnorm' estimates quantiles", "for a zero-modified lognormal distribution.  You have supplied an object", 
                "that assumes a different distribution."))
        class.x <- oldClass(x)
        if (!is.null(x$interval)) {
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        meanlog <- x$parameters["meanlog"]
        sdlog <- x$parameters["sdlog"]
        p.zero <- x$parameters["p.zero"]
        n <- x$sample.size
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
        if (n < 1 || any(x < 0)) 
            stop(paste("'x' must contain at least 1 non-missing value,", 
                "and all non-missing values of 'x' must be non-negative. ", 
                "This is not true for 'x' =", data.name))
        ret.list <- ezmlnorm(x, method = method)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        meanlog <- ret.list$parameters["meanlog"]
        sdlog <- ret.list$parameters["sdlog"]
        p.zero <- ret.list$parameters["p.zero"]
    }
    q <- qzmlnorm(p, meanlog = meanlog, sdlog = sdlog, p.zero = p.zero)
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
