ppccEvdGofTest <-
function (x, est.arg.list = NULL) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    n <- length(x)
    if (n < 10 || length(unique(x)) < 2) 
        stop(paste("'x' must contain at least 10 non-missing values,", 
            "and at least 2 distinct values. ", "This is not true for 'x' =", 
            data.name))
    if (n > 10000) 
        stop(paste("Too many observations.  Critical values of r", 
            "are available only for sample sizes between 10 and 10,000"))
    if (!is.null(est.arg.list) && !is.list(est.arg.list)) 
        stop("'est.arg.list' must be a list")
    ret.list <- do.call("eevd", c(list(x = x), est.arg.list))
    nrl <- names(ret.list)
    names(ret.list)[match("parameters", nrl)] <- "distribution.parameters"
    names(ret.list)[match("method", nrl)] <- "estimation.method"
    ret.list$data.name <- data.name
    ret.list$bad.obs <- bad.obs
    ret.list$dist.abb <- "evd"
    r <- cor(sort(x), qevd(ppoints(n, a = 0.44)))
    p <- ppccEvdGofTestPValue(r, n)
    ret.list <- c(ret.list, list(statistic = r, parameters = n, 
        p.value = p, alternative = paste("True cdf does not equal the\n", 
            space(33), "Extreme Value Distribution.", sep = ""), 
        method = "PPCC GOF", data = x))
    names(ret.list$statistic) <- "r"
    names(ret.list$parameters) <- "n"
    ret.list <- ret.list[c("distribution", "dist.abb", "distribution.parameters", 
        "n.param.est", "estimation.method", "statistic", "sample.size", 
        "parameters", "p.value", "alternative", "method", "data", 
        "data.name", "bad.obs")]
    oldClass(ret.list) <- "gof"
    ret.list
}
