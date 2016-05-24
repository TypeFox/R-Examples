egammaCensored <-
function (x, censored, method = "mle", censoring.side = "left", 
    ci = FALSE, ci.method = "profile.likelihood", ci.type = "two-sided", 
    conf.level = 0.95, n.bootstraps = 1000, pivot.statistic = "z", 
    ci.sample.size = sum(!censored)) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    if (!is.vector(censored, mode = "numeric") && !is.vector(censored, 
        mode = "logical")) 
        stop("'censored' must be a logical or numeric vector")
    if (length(censored) != length(x)) 
        stop("'censored' must be the same length as 'x'")
    data.name <- deparse(substitute(x))
    censoring.name <- deparse(substitute(censored))
    if ((bad.obs <- sum(!(ok <- is.finite(x) & is.finite(as.numeric(censored))))) > 
        0) {
        x <- x[ok]
        censored <- censored[ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' and 'censored' removed."))
    }
    if (is.numeric(censored)) {
        if (!all(censored == 0 | censored == 1)) 
            stop(paste("When 'censored' is a numeric vector, all values of", 
                "'censored' must be 0 (not censored) or 1 (censored)."))
        censored <- as.logical(censored)
    }
    if (any(x <= 0)) 
        stop("All values of 'x' (including censored ones) must be positive")
    n.cen <- sum(censored)
    if (n.cen == 0) 
        stop("No censored values indicated by 'censored'.")
    x.no.cen <- x[!censored]
    if (length(unique(x.no.cen)) < 2) 
        stop("'x' must contain at least 2 non-missing, uncensored, distinct values.")
    N <- length(x)
    method <- match.arg(method, c("mle"))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    x.cen <- x[censored]
    c.vec <- table(x.cen)
    cen.levels <- sort(unique(x.cen))
    ci.method <- match.arg(ci.method, c("normal.approx", "bootstrap", 
        "profile.likelihood"))
    ci.type <- match.arg(ci.type, c("two-sided", "lower", "upper"))
    if (ci && ci.method == "profile.likelihood" && method != 
        "mle") 
        stop("When ci.method=\"profile.likelihood\" you must set method=\"mle\"")
    pivot.statistic <- match.arg(pivot.statistic, c("z", "t"))
    if (conf.level <= 0 || conf.level >= 1) 
        stop("The value of 'conf.level' must be between 0 and 1.")
    est.fcn <- paste("egammaCensored", method, sep = ".")
    args.list <- switch(method, mle = list(x = x, censored = censored, 
        censoring.side = censoring.side, ci = ci, ci.method = ci.method, 
        ci.type = ci.type, conf.level = conf.level, ci.sample.size = ci.sample.size, 
        pivot.statistic = pivot.statistic))
    if (!ci || ci.method != "bootstrap") {
        param.ci.list <- do.call(est.fcn, args = args.list)
    }
    else {
        args.list$ci <- FALSE
        param.ci.list <- do.call(est.fcn, args = args.list)
        est.fcn <- paste("egammaAltCensored", method, sep = ".")
        param.list <- do.call(est.fcn, args = args.list)
        ci.list <- egammaAltCensored.bootstrap.ci(x = x, censored = censored, 
            censoring.side = censoring.side, est.fcn = est.fcn, 
            ci.type = ci.type, conf.level = conf.level, n.bootstraps = n.bootstraps, 
            obs.mean = param.list$parameters["mean"])
        param.ci.list <- c(param.ci.list, list(ci.obj = ci.list))
    }
    method.string <- switch(method, mle = "MLE")
    ret.list <- list(distribution = "Gamma", sample.size = N, 
        censoring.side = censoring.side, censoring.levels = cen.levels, 
        percent.censored = (100 * n.cen)/N, parameters = param.ci.list$parameters, 
        n.param.est = 2, method = method.string)
    ret.list <- c(ret.list, list(data.name = data.name, censoring.name = censoring.name, 
        bad.obs = bad.obs))
    if (ci) {
        ret.list <- c(ret.list, list(interval = param.ci.list$ci.obj))
        if (!is.null(param.ci.list$var.cov.params)) 
            ret.list <- c(ret.list, list(var.cov.params = param.ci.list$var.cov.params))
    }
    oldClass(ret.list) <- "estimateCensored"
    ret.list
}
