eexp <-
function (x, method = "mle/mme", ci = FALSE, ci.type = "two-sided", 
    ci.method = "exact", conf.level = 0.95) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    method <- match.arg(method)
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    n <- length(x)
    if (n < 1 || any(x < 0) || all(x == 0)) 
        stop(paste("'x' must contain at least one non-missing value,", 
            "all non-missing values of 'x' must be non-negative,", 
            "and at least one value of 'x' must be positive."))
    rate <- 1/mean(x)
    ret.list <- list(distribution = "Exponential", sample.size = n, 
        parameters = c(rate = rate), n.param.est = 1, method = method, 
        data.name = data.name, bad.obs = bad.obs)
    if (ci) {
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        ci.method <- match.arg(ci.method)
        if (conf.level <= 0 || conf.level >= 1) 
            stop("The value of 'conf.level' must be between 0 and 1.")
        ci.obj <- switch(ci.method, exact = ci.exp.exact(rate = rate, 
            n = n, ci.type = ci.type, alpha = 1 - conf.level))
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    oldClass(ret.list) <- "estimate"
    ret.list
}
