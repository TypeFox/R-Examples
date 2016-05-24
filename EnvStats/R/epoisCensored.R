epoisCensored <-
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
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    ci.type <- match.arg(ci.type, c("two-sided", "lower", "upper"))
    pivot.statistic <- match.arg(pivot.statistic, c("z", "t"))
    multiple <- TRUE
    T.vec <- unique(x[censored])
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
    args.list <- list(x = x, censored = censored, method = method, 
        censoring.side = censoring.side, ci = ci, ci.method = ci.method, 
        ci.type = ci.type, conf.level = conf.level, n.bootstraps = n.bootstraps, 
        ci.sample.size = ci.sample.size, pivot.statistic = pivot.statistic)
    if (multiple) {
        ret.list <- do.call("epoisMultiplyCensored", args = args.list)
    }
    else {
        ret.list <- do.call("epoisSinglyCensored", args = args.list)
    }
    ret.list$data.name <- data.name
    ret.list$censoring.name <- censoring.name
    ret.list
}
