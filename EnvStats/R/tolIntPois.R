tolIntPois <-
function (x, coverage = 0.95, cov.type = "content", ti.type = "two-sided", 
    conf.level = 0.95) 
{
    if (any(is.na(coverage))) 
        stop("Missing values not allowed in 'coverage'.")
    if (!is.numeric(coverage) || length(coverage) > 1 || coverage <= 
        0 || coverage >= 1) 
        stop("'coverage' must be a scalar greater than 0 and less than 1.")
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    ti.type <- match.arg(ti.type, c("two-sided", "lower", "upper"))
    if (!is.numeric(conf.level) || length(conf.level) > 1 || 
        conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be a scalar greater than 0 and less than 1.")
    if (x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored") {
        if (cov.type == "content") 
            stop(paste("When cov.type='content' you must supply", 
                "the original observations; you cannot supply the", 
                "result of calling 'epois'"))
        if (x$distribution != "Poisson") 
            stop(paste("'tolIntPois' creates tolerance intervals", 
                "for a Poisson distribution.  You have supplied an object", 
                "that assumes a different distribution."))
        class.x <- oldClass(x)
        if (!is.null(x$interval)) {
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        ret.list <- x
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
        n <- length(x)
        if (n < 2 || length(unique(x)) < 2) 
            stop("'x' must contain at least 2 finite, non-missing distinct values.")
        if (any(x < 0) || any(x != trunc(x))) 
            stop("All finite, non-missing values of 'x' must be non-negative integers")
        if (cov.type == "content") 
            ret.list <- epois(x, ci = TRUE, ci.method = "exact", 
                ci.type = ti.type, conf.level = conf.level)
        else ret.list <- epois(x)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
    }
    if (cov.type == "content") {
        conf.limits <- ret.list$interval$limits
        ret.list <- ret.list[-match("interval", names(ret.list))]
        switch(ti.type, `two-sided` = {
            ltl <- qpois((1 - coverage)/2, conf.limits["LCL"])
            utl <- qpois((1 + coverage)/2, conf.limits["UCL"])
        }, lower = {
            ltl <- qpois(1 - coverage, conf.limits["LCL"])
            ucl <- Inf
        }, upper = {
            ltl <- 0
            utl <- qpois(coverage, conf.limits["UCL"])
        })
        limits <- c(ltl, utl)
        names(limits) <- c("LTL", "UTL")
        ti.obj <- list(name = "Tolerance", coverage = coverage, 
            coverage.type = cov.type, limits = limits, type = ti.type, 
            method = "Zacks", conf.level = conf.level, sample.size = n)
    }
    else {
        ti.obj <- predIntPois(x, k = 1, pi.type = ti.type, conf.level = coverage)$interval
        ti.obj <- ti.obj[c("name", "limits", "type", "method", 
            "sample.size")]
        ti.obj$name <- "Tolerance"
        ti.obj$coverage <- coverage
        ti.obj$coverage.type <- cov.type
        names(ti.obj$limits) <- c("LTL", "UTL")
    }
    oldClass(ti.obj) <- "intervalEstimate"
    ret.list <- c(ret.list, list(interval = ti.obj))
    if (x.is.est.obj) 
        oldClass(ret.list) <- class.x
    else oldClass(ret.list) <- "estimate"
    ret.list
}
