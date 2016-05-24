eqlnorm <-
function (x, p = 0.5, method = "qmle", ci = FALSE, ci.method = "exact", 
    ci.type = "two-sided", conf.level = 0.95, digits = 0) 
{
    if (!is.vector(p, mode = "numeric") || is.factor(p)) 
        stop("'p' must be a numeric vector.")
    if (any(!is.finite(p))) 
        stop("NA/NaN/Inf values not allowed for 'p'.")
    if (any(p < 0) || any(p > 1)) 
        stop("All values of 'p' must be between 0 and 1.")
    method <- match.arg(method, c("qmle", "mvue"))
    ci.method <- match.arg(ci.method, c("exact", "normal.approx"))
    ci.type <- match.arg(ci.type, c("two-sided", "lower", "upper"))
    if (x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored") {
        if (x$distribution != "Lognormal") 
            stop(paste("'eqlnorm' estimates quantiles", "for a lognormal distribution.  You have supplied an object", 
                "that assumes a different distribution."))
        if (names(x$parameters[1]) == "mean") 
            stop(paste("You have suppled an object resulting from a call", 
                "to a function whose name begins with 'elnormAlt',", 
                "not 'elnorm'."))
        if (!is.null(x$interval)) {
            class.x <- oldClass(x)
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        new.x <- x
        names(new.x$parameters) <- c("mean", "sd")
        new.x$distribution <- "Normal"
        ret.obj <- eqnorm(new.x, p = p, ci = ci, ci.method = ci.method, 
            ci.type = ci.type, conf.level = conf.level, digits = digits)
        ret.obj$parameters <- x$parameters
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
        if (any(x <= 0)) 
            stop("All non-missing values of 'x' must be positive.")
        n <- length(x)
        if (n < 2 || length(unique(x)) < 2) 
            stop("'x' must contain at least 2 non-missing distinct values")
        ret.obj <- eqnorm(log(x), p = p, ci = ci, ci.method = ci.method, 
            ci.type = ci.type, conf.level = conf.level, digits = digits)
        ret.obj$data.name <- data.name
        ret.obj$bad.obs <- bad.obs
        names(ret.obj$parameters) <- c("meanlog", "sdlog")
    }
    ret.obj$distribution <- "Lognormal"
    ret.obj$quantiles <- exp(ret.obj$quantiles)
    if (ci) {
        ret.obj$interval$limits <- exp(ret.obj$interval$limits)
    }
    if (method == "mvue") {
        if (length(p) != 1 || p != 0.5) 
            stop(paste("The 'mvue' method is only available for", 
                "the median (i.e., p=0.5)"))
        meanlog <- ret.obj$parameters["meanlog"]
        sdlog <- ret.obj$parameters["sdlog"]
        s2 <- sdlog^2
        n <- ret.obj$sample.size
        df <- n - 1
        mhat <- exp(meanlog) * finneys.g(df, -s2/(2 * df))
        ret.obj$quantiles <- mhat
        names(ret.obj$quantiles) <- "Median"
        if (x.is.est.obj && x$method != "mvue") 
            ret.obj$quantile.method <- paste("quasi-mvue based on\n", 
                space(33), ret.obj$method, " Estimators", sep = "")
        else ret.obj$quantile.method <- "mvue"
        if (ci && ci.method == "normal.approx") {
            sd.mhat <- sqrt(exp(2 * meanlog) * ((finneys.g(df, 
                -s2/(2 * df))^2) - finneys.g(df, (-2 * s2)/df)))
            ret.obj$interval$limits <- ci.normal.approx(mhat, 
                sd.mhat, n, df, ci.type, alpha = 1 - conf.level)$limits
            ret.obj$interval$method <- "Normal Approx"
        }
    }
    ret.obj
}
