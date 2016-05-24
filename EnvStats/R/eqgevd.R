eqgevd <-
function (x, p = 0.5, method = "mle", pwme.method = "unbiased", 
    tsoe.method = "med", plot.pos.cons = c(a = 0.35, b = 0), 
    digits = 0) 
{
    if (!is.vector(p, mode = "numeric") || is.factor(p)) 
        stop("'p' must be a numeric vector.")
    if (any(!is.finite(p))) 
        stop("NA/NaN/Inf values not allowed in 'p'.")
    if (any(p < 0) || any(p > 1)) 
        stop("All values of 'p' must be between 0 and 1.")
    method <- match.arg(method, c("mle", "pwme", "tsoe"))
    pwme.method <- match.arg(pwme.method, c("unbiased", "plotting.position"))
    tsoe.method <- match.arg(tsoe.method, c("med", "lms", "lts"))
    if (x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored") {
        if (x$distribution != "Generalized Extreme Value") 
            stop(paste("'eqgevd' estimates quantiles", "for a generalized extreme value distribution.  You have supplied an object", 
                "that assumes a different distribution."))
        class.x <- oldClass(x)
        if (!is.null(x$interval)) {
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        location <- x$parameters["location"]
        scale <- x$parameters["scale"]
        shape <- x$parameters["shape"]
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
        if (n < 3 || length(unique(x)) < 3) 
            stop(paste("'x' must contain at least 3 non-missing distinct values. ", 
                "This is not true for 'x' =", data.name))
        ret.list <- egevd(x, method = method, pwme.method = pwme.method, 
            tsoe.method = tsoe.method, plot.pos.cons = plot.pos.cons)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        location <- ret.list$parameters["location"]
        scale <- ret.list$parameters["scale"]
        shape <- ret.list$parameters["shape"]
    }
    q <- qgevd(p, location = location, scale = scale, shape = shape)
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
