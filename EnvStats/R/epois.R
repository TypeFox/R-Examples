epois <-
function (x, method = "mle/mme/mvue", ci = FALSE, ci.type = "two-sided", 
    ci.method = "exact", conf.level = 0.95) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    method <- match.arg(method)
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    n <- length(x)
    if (n < 1 || any(x < 0) || any(x != trunc(x))) 
        stop(paste("x must contain at least one non-missing value,", 
            "and all non-missing values of x must be non-negative integers. ", 
            "This is not true for x =", data.name))
    dist.param <- c(lambda = mean(x))
    ret.list <- list(distribution = "Poisson", sample.size = n, 
        parameters = dist.param, n.param.est = 1, method = method, 
        data.name = data.name, bad.obs = bad.obs)
    if (ci) {
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        ci.method <- match.arg(ci.method, c("exact", "pearson.hartley.approx", 
            "normal.approx"))
        if (conf.level <= 0 || conf.level >= 1) 
            stop("The value of 'conf.level' must be between 0 and 1.")
        if (dist.param == 0 && ci.method != "exact") 
            stop(paste("All non-missing values in 'x' are 0. ", 
                "You must use ci.method='exact' in this case."))
        ci.obj <- switch(ci.method, exact = ci.pois.exact(x, 
            alpha = 1 - conf.level, ci.type = ci.type), pearson.hartley.approx = ci.pois.pearson.hartley.approx(x, 
            alpha = 1 - conf.level, ci.type = ci.type), normal.approx = ci.normal.approx(theta.hat = dist.param, 
            sd.theta.hat = sqrt(dist.param/n), n = n, df = Inf, 
            ci.type = ci.type, alpha = 1 - conf.level, lb = 0))
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    oldClass(ret.list) <- "estimate"
    ret.list
}
