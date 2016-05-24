eqlnormCensored <-
function (x, censored, censoring.side = "left", p = 0.5, method = "mle", 
    ci = FALSE, ci.method = "exact.for.complete", ci.type = "two-sided", 
    conf.level = 0.95, digits = 0, nmc = 1000, seed = NULL) 
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
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    if (!is.vector(p, mode = "numeric")) 
        stop("'p' must be a numeric vector.")
    if (any(!is.finite(p))) 
        stop("NA/NaN/Inf values not allowed in 'p'.")
    if (any(p < 0) || any(p > 1)) 
        stop("All values of 'p' must be between 0 and 1")
    ci.method <- match.arg(ci.method, c("exact.for.complete", 
        "gpq", "normal.approx"))
    ret.list <- eqnormCensored(x = log(x), censored = censored, 
        censoring.side = censoring.side, p = p, method = method, 
        ci = ci, ci.method = ci.method, ci.type = ci.type, conf.level = conf.level, 
        digits = digits, seed = seed, nmc = nmc)
    ret.list$data.name <- data.name
    ret.list$censoring.name <- censoring.name
    ret.list$bad.obs <- bad.obs
    ret.list$distribution <- "Lognormal"
    ret.list$censoring.levels <- exp(ret.list$censoring.levels)
    names(ret.list$parameters) <- c("meanlog", "sdlog")
    ret.list$quantiles <- exp(ret.list$quantiles)
    if (ci) 
        ret.list$interval$limits <- exp(ret.list$interval$limits)
    ret.list
}
