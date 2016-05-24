enorm <-
function (x, method = "mvue", ci = FALSE, ci.type = "two-sided", 
    ci.method = "exact", conf.level = 0.95, ci.param = "mean") 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if (length(ci) != 1 || !is.logical(ci)) 
        stop("The argument 'ci' must be a logical scalar")
    if (ci) 
        method <- "mvue"
    else method <- match.arg(method, c("mvue", "mle/mme"))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    n <- length(x)
    if (n < 2 || length(unique(x)) < 2) 
        stop("'x' must contain at least 2 non-missing distinct values")
    muhat <- mean(x)
    sdhat <- ifelse(method == "mvue", sd(x), sqrt((n - 1)/n) * 
        sd(x))
    ret.list <- list(distribution = "Normal", sample.size = n, 
        parameters = c(mean = muhat, sd = sdhat), n.param.est = 2, 
        method = method, data.name = data.name, bad.obs = bad.obs)
    if (ci) {
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        ci.method <- match.arg(ci.method, "exact")
        if (conf.level <= 0 || conf.level >= 1) 
            stop("The value of 'conf.level' must be between 0 and 1.")
        ci.param <- match.arg(ci.param, c("mean", "variance"))
        ci.obj <- switch(ci.param, mean = {
            switch(ci.method, exact = ci.norm(muhat, sdhat, n, 
                ci.type, alpha = 1 - conf.level))
        }, variance = {
            switch(ci.method, exact = ci.norm.var(sdhat = sdhat, 
                n = n, ci.type = ci.type, alpha = 1 - conf.level))
        })
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    oldClass(ret.list) <- "estimate"
    ret.list
}
