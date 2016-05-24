elnormSinglyCensored <-
function (x, censored, method = "mle", censoring.side = "left", 
    ci = FALSE, ci.method = "profile.likelihood", ci.type = "two-sided", 
    conf.level = 0.95, n.bootstraps = 1000, pivot.statistic = "z", 
    nmc = 1000, seed = NULL, ...) 
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
    n.cen <- sum(censored)
    if (n.cen == 0) 
        stop("No censored values indicated by 'censored'.")
    if (any(x <= 0)) 
        stop("All values of 'x' (including censored ones) must be positive")
    x.no.cen <- x[!censored]
    if (length(unique(x.no.cen)) < 2) 
        stop("'x' must contain at least 2 non-missing, uncensored, distinct values.")
    N <- length(x)
    method <- match.arg(method, c("mle", "bcmle", "qq.reg", "qq.reg.w.cen.level", 
        "impute.w.qq.reg", "impute.w.qq.reg.w.cen.level", "impute.w.mle", 
        "iterative.impute.w.qq.reg", "half.cen.level", "m.est"))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    T1 <- unique(x[censored])
    if (length(T1) > 1) 
        stop(paste("More than one censoring level.  Use 'elnormMultiplyCensored'."))
    if (censoring.side == "left") {
        if (T1 > min(x.no.cen)) 
            stop(paste("For singly left-censored data,", "all uncensored observations must be bigger than", 
                "or equal to the censoring level. ", "Use elnormMultiplyCensored."))
    }
    else {
        if (T1 < max(x.no.cen)) 
            stop(paste("For singly right-censored data,", "all uncensored observations must be less than", 
                "or equal to the censoring level. ", "Use elnormMultiplyCensored."))
    }
    if (method == "half.cen.level" && censoring.side == "right") 
        stop(paste("The method 'half.cen.level' is applicable only for", 
            "left-censored data with a positive censoring level"))
    ci.method <- match.arg(ci.method, c("normal.approx", "normal.approx.w.cov", 
        "bootstrap", "profile.likelihood", "gpq"))
    ci.type <- match.arg(ci.type, c("two-sided", "lower", "upper"))
    if (ci && ci.method == "profile.likelihood" && method != 
        "mle") 
        stop("When ci.method=\"profile.likelihood\" you must set method=\"mle\"")
    pivot.statistic <- match.arg(pivot.statistic, c("z", "t"))
    ret.list <- enormSinglyCensored(x = log(x), censored = censored, 
        method = method, censoring.side = censoring.side, ci = ci, 
        ci.method = ci.method, ci.type = ci.type, conf.level = conf.level, 
        n.bootstraps = n.bootstraps, pivot.statistic = pivot.statistic, 
        nmc = nmc, seed = seed, ...)
    ret.list$distribution <- "Lognormal"
    ret.list$censoring.levels <- T1
    names(ret.list$parameters) <- c("meanlog", "sdlog")
    ret.list$data.name <- data.name
    ret.list$censoring.name <- censoring.name
    if (!is.null(ret.list$interval)) 
        ret.list$interval$parameter <- "meanlog"
    ret.list
}
