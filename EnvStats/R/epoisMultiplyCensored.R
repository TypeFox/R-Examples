epoisMultiplyCensored <-
function (x, censored, method = "mle", censoring.side = "left", 
    ci = FALSE, ci.method = "profile.likelihood", ci.type = "two-sided", 
    conf.level = 0.95, n.bootstraps = 1000, ci.sample.size = N - 
        n.cen, pivot.statistic = "z") 
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
    N <- length(x)
    if (N < 1) 
        stop("No valid values in 'x' and 'censored'")
    n.cen <- sum(censored)
    if (n.cen == 0) 
        stop("No censored values indicated by 'censored'.")
    x.no.cen <- x[!censored]
    if ((N - n.cen) < 1 || any(x < 0) || any(x != trunc(x))) 
        stop(paste("'x' must contain at least one non-missing,", 
            "uncensored value, and all non-missing values of 'x'", 
            "(including censored observations) must be non-negative integers."))
    method <- match.arg(method, c("mle", "half.cen.level"))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    x.cen <- x[censored]
    c.vec <- table(x.cen)
    cen.levels <- sort(unique(x.cen))
    K <- length(c.vec)
    if (method == "half.cen.level" && (censoring.side == "right" || 
        any(cen.levels < 2 * .Machine$double.eps))) 
        stop(paste("The method 'half.cen.level' is applicable only for", 
            "left-censored data with positive censoring levels"))
    ci.method <- match.arg(ci.method, c("normal.approx", "bootstrap", 
        "profile.likelihood"))
    ci.type <- match.arg(ci.type, c("two-sided", "lower", "upper"))
    if (ci & ci.method == "profile.likelihood" && method != "mle") 
        stop("When ci.method=\"profile.likelihood\" you must set method=\"mle\"")
    pivot.statistic <- match.arg(pivot.statistic, c("z", "t"))
    if (ci && ci.method == "normal.approx" && pivot.statistic == 
        "t" && (ci.sample.size < 2 | ci.sample.size > N)) 
        stop("'ci.sample.size' must be between 2 and the number of non-missing observations")
    est.fcn <- paste("epoisMultiplyCensored", method, sep = ".")
    args.list <- list(x = x, censored = censored, censoring.side = censoring.side, 
        ci = ci, ci.method = ci.method, ci.type = ci.type, conf.level = conf.level, 
        ci.sample.size = ci.sample.size, pivot.statistic = pivot.statistic)
    if (!ci || ci.method != "bootstrap") {
        param.ci.list <- do.call(est.fcn, args = args.list)
    }
    else {
        args.list$ci <- FALSE
        param.ci.list <- do.call(est.fcn, args = args.list)
        ci.list <- epoisMultiplyCensored.bootstrap.ci(x = x, 
            censored = censored, censoring.side = censoring.side, 
            est.fcn = est.fcn, ci.type = ci.type, conf.level = conf.level, 
            n.bootstraps = n.bootstraps, obs.lambda = param.ci.list$parameters["lambda"])
        param.ci.list <- c(param.ci.list, list(ci.obj = ci.list))
    }
    method.string <- switch(method, mle = "MLE", half.cen.level = "Half Censoring Level")
    ret.list <- list(distribution = "Poisson", sample.size = N, 
        censoring.side = censoring.side, censoring.levels = cen.levels, 
        percent.censored = (100 * n.cen)/N, parameters = param.ci.list$parameters, 
        n.param.est = 1, method = method.string, data.name = data.name, 
        censoring.name = censoring.name, bad.obs = bad.obs)
    if (ci) {
        ret.list <- c(ret.list, list(interval = param.ci.list$ci.obj))
        if (!is.null(param.ci.list$var.cov.params)) 
            ret.list <- c(ret.list, list(var.cov.params = param.ci.list$var.cov.params))
    }
    oldClass(ret.list) <- "estimateCensored"
    ret.list
}
