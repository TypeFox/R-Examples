predIntLnormSimultaneous <-
function (x, n.geomean = 1, k = 1, m = 2, r = 1, rule = "k.of.m", 
    delta.over.sigma = 0, pi.type = "upper", conf.level = 0.95, 
    K.tol = .Machine$double.eps^0.5) 
{
    rule <- match.arg(rule, c("k.of.m", "CA", "Modified.CA"))
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    switch(rule, k.of.m = {
        if (!is.vector(n.geomean, mode = "numeric") || length(n.geomean) != 
            1 || n.geomean != trunc(n.geomean) || n.geomean < 
            1 || !is.vector(k, mode = "numeric") || length(k) != 
            1 || k != trunc(k) || k < 1 || !is.vector(m, mode = "numeric") || 
            length(m) != 1 || m != trunc(m) || m < 1 || !is.vector(r, 
            mode = "numeric") || length(r) != 1 || r != trunc(r) || 
            r < 1 || k > m) stop(paste("'n.geomean', 'k', 'm', and 'r' must be positive integers,", 
            "and 'k' must be between 1 and 'm'"))
    }, CA = {
        if (!is.vector(n.geomean, mode = "numeric") || length(n.geomean) != 
            1 || n.geomean != trunc(n.geomean) || n.geomean < 
            1 || !is.vector(m, mode = "numeric") || length(m) != 
            1 || m != trunc(m) || m < 1 || !is.vector(r, mode = "numeric") || 
            length(r) != 1 || r != trunc(r) || r < 1) stop("'n.geomean', 'm', and 'r' must be positive integers")
    }, Modified.CA = {
        if (!is.vector(n.geomean, mode = "numeric") || length(n.geomean) != 
            1 || n.geomean != trunc(n.geomean) || n.geomean < 
            1 || !is.vector(m, mode = "numeric") || length(m) != 
            1 || m != trunc(m) || m < 1 || !is.vector(r, mode = "numeric") || 
            length(r) != 1 || r != trunc(r) || r < 1) stop("'n.geomean', 'm', and 'r' must be positive integers")
        m <- 4
    })
    if (!is.vector(delta.over.sigma, mode = "numeric") || length(delta.over.sigma) != 
        1 || !is.finite(delta.over.sigma)) 
        stop("'delta.over.sigma' must be a finite numeric scalar.")
    if (!is.vector(conf.level, mode = "numeric") || length(conf.level) != 
        1 || conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be a scalar greater than 0 and less than 1.")
    if (x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored") {
        if (x$distribution != "Lognormal") 
            stop(paste("'predIntLnormSimultaneous' creates prediction intervals", 
                "for a Lognormal distribution.  You have supplied an object", 
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
        ret.list <- predIntNormSimultaneous(new.x, n.mean = n.geomean, 
            k = k, m = m, r = r, rule = rule, delta.over.sigma = delta.over.sigma, 
            pi.type = pi.type, conf.level = conf.level, K.tol = K.tol)
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
        ret.list <- predIntNormSimultaneous(log(x), n.mean = n.geomean, 
            k = k, m = m, r = r, rule = rule, delta.over.sigma = delta.over.sigma, 
            pi.type = pi.type, conf.level = conf.level, K.tol = K.tol)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        names(ret.list$parameters) <- c("meanlog", "sdlog")
    }
    names(ret.list$interval)[names(ret.list$interval) == "n.mean"] <- "n.geomean"
    ret.list$distribution <- "Lognormal"
    ret.list$interval$limits <- exp(ret.list$interval$limits)
    ret.list
}
