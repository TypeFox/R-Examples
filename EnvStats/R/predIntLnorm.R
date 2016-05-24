predIntLnorm <-
function (x, n.geomean = 1, k = 1, method = "Bonferroni", pi.type = "two-sided", 
    conf.level = 0.95) 
{
    if (!is.vector(k, mode = "numeric") || length(k) != 1 || 
        !is.vector(n.geomean, mode = "numeric") || length(n.geomean) != 
        1 || !is.vector(conf.level, mode = "numeric") || length(conf.level) != 
        1) 
        stop("'k', 'n.geomean', and 'conf.level' must be numeric scalars")
    if (k != trunc(k) || k < 1) 
        stop("'k' must be an integer greater than 0")
    if (n.geomean != trunc(n.geomean) || n.geomean < 1) 
        stop("'n.geomean' must be an integer greater than 0")
    if (conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be greater than 0 and less than 1.")
    method <- match.arg(method, c("Bonferroni", "exact"))
    pi.type <- match.arg(pi.type, c("two-sided", "lower", "upper"))
    if (data.class(x) == "estimate" || data.class(x) == "estimateCensored") {
        if (x$distribution != "Lognormal") 
            stop(paste("'predIntLnorm' creates prediction intervals", 
                "for a lognormal distribution.  You have supplied an object", 
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
        ret.list <- predIntNorm(new.x, n.mean = n.geomean, k = k, 
            method = method, pi.type = pi.type, conf.level = conf.level)
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
        if (any(x <= 0)) 
            stop("All non-missing values of 'x' must be positive")
        n <- length(x)
        if (n < 2 || length(unique(x)) < 2) 
            stop(paste("'x' must contain at least 2 non-missing distinct values. ", 
                "This is not true for 'x' =", data.name))
        ret.list <- predIntNorm(log(x), n.mean = n.geomean, k = k, 
            method = method, pi.type = pi.type, conf.level = conf.level)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        names(ret.list$parameters) <- c("meanlog", "sdlog")
    }
    names(ret.list$interval)[names(ret.list$interval) == "n.mean"] <- "n.geomean"
    ret.list$distribution <- "Lognormal"
    ret.list$interval$limits <- exp(ret.list$interval$limits)
    ret.list
}
