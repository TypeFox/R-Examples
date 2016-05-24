eqnorm <-
function (x, p = 0.5, method = "qmle", ci = FALSE, ci.method = "exact", 
    ci.type = "two-sided", conf.level = 0.95, digits = 0, warn = TRUE) 
{
    if (!is.vector(p, mode = "numeric") || is.factor(p)) 
        stop("'p' must be a numeric vector.")
    if (any(!is.finite(p))) 
        stop("NA/NaN/Inf values not allowed in 'p'.")
    if (any(p < 0) || any(p > 1)) 
        stop("All values of 'p' must be between 0 and 1.")
    method <- match.arg(method, c("qmle"))
    ci.method <- match.arg(ci.method, c("exact", "normal.approx"))
    if (x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored") {
        if (x$distribution != "Normal") 
            stop(paste("'eqnorm' estimates quantiles", "for a normal distribution.  You have supplied an object", 
                "that assumes a different distribution."))
        class.x <- oldClass(x)
        if (!is.null(x$interval)) {
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        xbar <- x$parameters["mean"]
        s <- x$parameters["sd"]
        n <- x$sample.size
        ret.list <- x
        if (ci && ci.method == "exact" && ret.list$method != 
            "mvue" && warn) 
            warning(paste("When ci=T and ci.method=\"exact\", the supplied object", 
                "'x' that is of class 'estimate' should have used", 
                "method=\"mvue\" for estimation.\n"))
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
        ret.list <- enorm(x, method = "mvue")
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        xbar <- ret.list$parameters["mean"]
        s <- ret.list$parameters["sd"]
    }
    q <- qnorm(p, mean = xbar, sd = s)
    if (length(p) == 1 && p == 0.5) 
        names(q) <- "Median"
    else {
        pct <- round(100 * p, digits)
        names(q) <- paste(pct, number.suffix(pct), " %ile", sep = "")
    }
    ret.list <- c(ret.list, list(quantiles = q))
    if (x.is.est.obj && x$method != "mvue") 
        ret.list$quantile.method <- paste("Quantile(s) Based on\n", 
            space(33), ret.list$method, " Estimators", sep = "")
    else ret.list$quantile.method <- "qmle"
    if (ci) {
        if (length(p) > 1 || p <= 0 || p >= 1) 
            stop(paste("When 'ci' = TRUE, 'p' must be a scalar", 
                "larger than 0 and less than 1."))
        if (p == 0.5) 
            ci.method <- "exact"
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        if (conf.level <= 0 || conf.level >= 1) 
            stop("The value of 'conf.level' must be between 0 and 1.")
        ci.obj <- ci.qnorm(p = p, muhat = xbar, sdhat = s, n = n, 
            method = ci.method, ci.type = ci.type, alpha = 1 - 
                conf.level, digits = digits)
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    if (x.is.est.obj) 
        oldClass(ret.list) <- class.x
    else oldClass(ret.list) <- "estimate"
    ret.list
}
