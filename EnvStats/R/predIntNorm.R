predIntNorm <-
function (x, n.mean = 1, k = 1, method = "Bonferroni", pi.type = "two-sided", 
    conf.level = 0.95) 
{
    if (!is.vector(k, mode = "numeric") || length(k) != 1 || 
        k != trunc(k) || k < 1 || !is.vector(n.mean, mode = "numeric") || 
        length(n.mean) != 1 || n.mean != trunc(n.mean) || n.mean < 
        1) 
        stop("'k' and 'n.mean' must be positive integers")
    method <- match.arg(method, c("Bonferroni", "exact"))
    pi.type <- match.arg(pi.type, c("two-sided", "lower", "upper"))
    if (!is.vector(conf.level, mode = "numeric") || length(conf.level) != 
        1 || conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be a scalar greater than 0 and less than 1.")
    if (x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored") {
        if (x$distribution != "Normal") 
            stop(paste("'predIntNorm' creates prediction intervals", 
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
    if (k == 1) 
        method <- "exact"
    df <- n - 1
    K <- predIntNormK(n = n, df = df, n.mean = n.mean, k = k, 
        method = method, pi.type = pi.type, conf.level = conf.level)
    switch(pi.type, `two-sided` = {
        LPL <- xbar - K * s
        UPL <- xbar + K * s
    }, lower = {
        LPL <- xbar - K * s
        UPL <- Inf
    }, upper = {
        LPL <- -Inf
        UPL <- xbar + K * s
    })
    limits <- c(LPL, UPL)
    names(limits) <- c("LPL", "UPL")
    pi.obj <- list(name = "Prediction", limits = limits, type = pi.type, 
        method = method, conf.level = conf.level, sample.size = n, 
        dof = df, n.mean = n.mean, k = k)
    oldClass(pi.obj) <- "intervalEstimate"
    ret.list <- c(ret.list, list(interval = pi.obj))
    if (x.is.est.obj) 
        oldClass(ret.list) <- class.x
    else oldClass(ret.list) <- "estimate"
    ret.list
}
