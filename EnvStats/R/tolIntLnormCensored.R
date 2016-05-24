tolIntLnormCensored <-
function (x, censored, censoring.side = "left", coverage = 0.95, 
    cov.type = "content", ti.type = "two-sided", conf.level = 0.95, 
    method = "mle", ti.method = "exact.for.complete", seed = NULL, 
    nmc = 1000) 
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
        stop("'x' must contain at least 2 non-missing, uncensored, distinct values")
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    if (!is.numeric(coverage) || length(coverage) != 1 || is.na(coverage) || 
        coverage <= 0 || coverage >= 1) 
        stop("'coverage' must be a scalar greater than 0 and less than 1.")
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    ti.type <- match.arg(ti.type, c("two-sided", "lower", "upper"))
    if (!is.numeric(conf.level) || length(conf.level) > 1 || 
        conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be a scalar greater than 0 and less than 1.")
    multiple <- TRUE
    T.vec <- sort(unique(x[censored]))
    if (length(T.vec) == 1) {
        if (censoring.side == "left") {
            if (T.vec <= min(x.no.cen)) 
                multiple <- FALSE
        }
        else {
            if (T.vec >= max(x.no.cen)) 
                multiple <- FALSE
        }
    }
    if (multiple) {
        method <- match.arg(method, c("mle", "qq.reg", "impute.w.qq.reg", 
            "half.cen.level"))
    }
    else {
        method <- match.arg(method, c("mle", "bcmle", "qq.reg", 
            "qq.reg.w.cen.level", "impute.w.qq.reg", "impute.w.qq.reg.w.cen.level", 
            "impute.w.mle", "iterative.impute.w.qq.reg", "m.est", 
            "half.cen.level"))
    }
    ti.method <- match.arg(ti.method, c("exact.for.complete", 
        "wald.wolfowitz.for.complete", "gpq"))
    if (ti.method == "gpq" && cov.type != "content") 
        stop("When ti.method='gpq' you must set cov.type='content'")
    ret.list <- tolIntNormCensored(x = log(x), censored = censored, 
        censoring.side = censoring.side, coverage = coverage, 
        cov.type = cov.type, ti.type = ti.type, conf.level = conf.level, 
        method = method, ti.method = ti.method, seed = seed, 
        nmc = nmc)
    names(ret.list$parameters) <- c("meanlog", "sdlog")
    ret.list$censoring.levels <- T.vec
    ret.list$data.name <- data.name
    ret.list$bad.obs <- bad.obs
    ret.list$distribution <- "Lognormal"
    ret.list$interval$limits <- exp(ret.list$interval$limits)
    ret.list
}
