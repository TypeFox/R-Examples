elnormMultiplyCensored <-
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
    if (any(x <= 0)) 
        stop("All values of 'x' (including censored ones) must be positive")
    n.cen <- sum(censored)
    if (n.cen == 0) 
        stop("No censored values indicated by 'censored'.")
    x.no.cen <- x[!censored]
    if (length(unique(x.no.cen)) < 2) 
        stop("'x' must contain at least 2 non-missing, uncensored, distinct values.")
    N <- length(x)
    method <- match.arg(method, c("mle", "qq.reg", "impute.w.qq.reg", 
        "half.cen.level"))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    if (method == "half.cen.level" && censoring.side == "right") 
        stop(paste("The method 'half.cen.level' is applicable only for", 
            "left-censored data with positive censoring levels"))
    ci.method <- match.arg(ci.method, c("normal.approx", "bootstrap", 
        "profile.likelihood", "gpq"))
    ci.type <- match.arg(ci.type, c("two-sided", "lower", "upper"))
    if (ci && ci.method == "profile.likelihood" && method != 
        "mle") 
        stop("When ci.method=\"profile.likelihood\" you must set method=\"mle\"")
    pivot.statistic <- match.arg(pivot.statistic, c("z", "t"))
    x.cen <- x[censored]
    cen.levels <- sort(unique(x.cen))
    ret.list <- enormMultiplyCensored(x = log(x), censored = censored, 
        method = method, censoring.side = censoring.side, ci = ci, 
        ci.method = ci.method, ci.type = ci.type, conf.level = conf.level, 
        n.bootstraps = n.bootstraps, pivot.statistic = pivot.statistic, 
        nmc = nmc, seed = seed, ...)
    ret.list$distribution <- "Lognormal"
    ret.list$censoring.levels <- cen.levels
    names(ret.list$parameters) <- c("meanlog", "sdlog")
    ret.list$data.name <- data.name
    ret.list$censoring.name <- censoring.name
    if (!is.null(ret.list$interval)) 
        ret.list$interval$parameter <- "meanlog"
    ret.list
}
