tolIntNorm <-
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
    if (x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored") {
        if (x$distribution != "Normal") 
            stop(paste("'tolIntNorm' creates tolerance intervals", 
                "for a normal distribution.  You have supplied an object", 
                "that assumes a different distribution."))
        class.x <- oldClass(x)
        if (!is.null(x$interval)) {
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        xbar <- x$parameters["mean"]
        s <- x$parameters["sd"]
        n <- x$sample.size
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
            stop(paste("'x' must contain at least 2 non-missing distinct values. ", 
                "This is not true for 'x' =", data.name))
        ret.list <- enorm(x)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        xbar <- ret.list$parameters["mean"]
        s <- ret.list$parameters["sd"]
    }
    df <- n - 1
    K <- tolIntNormK(n = n, df = df, coverage = coverage, cov.type = cov.type, 
        ti.type = ti.type, conf.level = conf.level, method = method)
    hw <- K * s
    switch(ti.type, `two-sided` = {
        ltl <- xbar - hw
        utl <- xbar + hw
    }, lower = {
        ltl <- xbar - hw
        utl <- Inf
    }, upper = {
        ltl <- -Inf
        utl <- xbar + hw
    })
    method <- ifelse(cov.type == "content" && ti.type == "two-sided" && 
        method == "wald.wolfowitz", "Wald-Wolfowitz Approx", 
        "Exact")
    limits <- c(ltl, utl)
    names(limits) <- c("LTL", "UTL")
    if (cov.type == "content") 
        ti.obj <- list(name = "Tolerance", coverage = coverage, 
            coverage.type = cov.type, limits = limits, type = ti.type, 
            method = method, conf.level = conf.level, sample.size = n, 
            dof = df)
    else ti.obj <- list(name = "Tolerance", coverage = coverage, 
        coverage.type = cov.type, limits = limits, type = ti.type, 
        method = method, sample.size = n, dof = df)
    oldClass(ti.obj) <- "intervalEstimate"
    ret.list <- c(ret.list, list(interval = ti.obj))
    if (x.is.est.obj) 
        oldClass(ret.list) <- class.x
    else oldClass(ret.list) <- "estimate"
    ret.list
}
