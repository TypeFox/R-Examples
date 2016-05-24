tolIntLnormAlt <-
function (x, coverage = 0.95, cov.type = "content", ti.type = "two-sided", 
    conf.level = 0.95, method = "exact", est.method = "mvue") 
{
    if (any(is.na(coverage))) 
        stop("Missing values not allowed in 'coverage'.")
    if (!is.numeric(coverage) || length(coverage) > 1 || coverage <= 
        0 || coverage >= 1) 
        stop("'coverage' must be a scalar greater than 0 and less than 1.")
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    ti.type <- match.arg(ti.type, c("two-sided", "lower", "upper"))
    method <- match.arg(method, c("exact", "wald.wolfowitz"))
    if (!is.numeric(conf.level) || length(conf.level) > 1 || 
        conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be a scalar greater than 0 and less than 1.")
    if (!is.vector(x, mode = "numeric")) 
        stop(paste("'x' must be a numeric vector"))
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    if (any(x < 0)) 
        stop("All finite, non-missing values of 'x' must be positive")
    n <- length(x)
    if (n < 2 || length(unique(x)) < 2) 
        stop(paste("'x' must contain at least 2 non-missing distinct values. ", 
            "This is not true for 'x' =", data.name))
    ret.list <- tolIntNorm(log(x), coverage = coverage, cov.type = cov.type, 
        ti.type = ti.type, conf.level = conf.level, method = method)
    est.list <- do.call("elnormAlt", args = list(x = x, method = est.method))
    ret.list$parameters <- est.list$parameters
    ret.list$method <- est.list$method
    ret.list$data.name <- data.name
    ret.list$bad.obs <- bad.obs
    ret.list$distribution <- "Lognormal"
    ret.list$interval$limits <- exp(ret.list$interval$limits)
    ret.list
}
