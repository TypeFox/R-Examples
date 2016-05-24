tolIntLnorm <-
function (x, coverage = 0.95, cov.type = "content", ti.type = "two-sided", 
    conf.level = 0.95, method = "exact") 
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
    if (data.class(x) == "estimate" || data.class(x) == "estimateCensored") {
        if (x$distribution != "Lognormal") 
            stop(paste("'tolIntLnorm' creates tolerance intervals", 
                "for a normal distribution.  You have supplied an object", 
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
        ret.list <- tolIntNorm(new.x, coverage = coverage, cov.type = cov.type, 
            ti.type = ti.type, conf.level = conf.level, method = method)
        ret.list$parameters <- x$parameters
    }
    else {
        if (!is.vector(x, mode = "numeric")) 
            stop(paste("'x' must be either a list that inherits from", 
                "the class 'estimate', or else a numeric vector"))
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
        names(ret.list$parameters) <- c("meanlog", "sdlog")
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
    }
    ret.list$distribution <- "Lognormal"
    ret.list$interval$limits <- exp(ret.list$interval$limits)
    ret.list
}
