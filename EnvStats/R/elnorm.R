elnorm <-
function (x, method = "mvue", ci = FALSE, ci.type = "two-sided", 
    ci.method = "exact", conf.level = 0.95) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector.")
    data.name <- deparse(substitute(x))
    method <- ifelse(ci, "mvue", match.arg(method, c("mvue", 
        "mle/mme")))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    n <- length(x)
    if (n < 2 || length(unique(x)) < 2 || any(x <= 0)) 
        stop(paste("'x' must contain at least 2 non-missing distinct values,", 
            "and all non-missing values must be positive."))
    log.x <- log(x)
    muhat <- mean(log.x)
    sdhat <- ifelse(method == "mvue", sd(log.x), sqrt((n - 1)/n) * 
        sd(log.x))
    ret.list <- list(distribution = "Lognormal", sample.size = n, 
        parameters = c(meanlog = muhat, sdlog = sdhat), n.param.est = 2, 
        method = method, data.name = data.name, bad.obs = bad.obs)
    if (ci) {
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        ci.method <- match.arg(ci.method, "exact")
        if (conf.level <= 0 || conf.level >= 1) 
            stop("The value of 'conf.level' must be between 0 and 1.")
        ci.obj <- switch(ci.method, exact = ci.norm(muhat, sdhat, 
            n, ci.type, alpha = 1 - conf.level))
        ci.obj$parameter <- "meanlog"
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    oldClass(ret.list) <- "estimate"
    ret.list
}
