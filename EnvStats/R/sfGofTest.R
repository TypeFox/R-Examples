sfGofTest <-
function (x, distribution = c("norm", "lnorm", "lnormAlt", "zmnorm", 
    "zmlnorm", "zmlnormAlt"), est.arg.list = NULL) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    distribution <- match.arg(distribution)
    if (any(distribution == c("lnorm", "lnormAlt")) && any(x <= 
        0)) 
        stop("All values of 'x' must be positive for a lognormal distribution")
    if (any(distribution == c("zmlnorm", "zmlnormAlt")) && any(x < 
        0)) 
        stop(paste("All values of 'x' must be non-negative for a", 
            "zero-modified lognormal distribution"))
    n <- switch(distribution, zmnorm = sum(x != 0), zmlnorm = , 
        zmlnormAlt = sum(x > 0), length(x))
    if (n < 5 || length(unique(x)) < 2) 
        stop(paste("'x' must contain at least 5 non-missing values,", 
            "and at least 2 distinct values. ", "This is not true for 'x' =", 
            data.name))
    if (n > 5000) 
        warning(paste("Too many observations.  This approximation only works", 
            "if the number of observations is between 5 and 5000"))
    est.fcn <- paste("e", distribution, sep = "")
    ret.list <- do.call(est.fcn, c(list(x = x), est.arg.list))
    nrl <- names(ret.list)
    names(ret.list)[match("parameters", nrl)] <- "distribution.parameters"
    names(ret.list)[match("method", nrl)] <- "estimation.method"
    ret.list$data.name <- data.name
    ret.list$bad.obs <- bad.obs
    ret.list$dist.abb <- distribution
    new.x <- switch(distribution, norm = sort(x), lnorm = , lnormAlt = sort(log(x)), 
        zmnorm = sort(x[x != 0]), zmlnorm = , zmlnormAlt = sort(log(x[x > 
            0])))
    W <- ppccNorm(new.x)^2
    w <- log(1 - W)
    nu <- log(n)
    y <- log(nu) - nu
    mu <- -1.2725 + 1.0521 * y
    y <- log(nu) + 2/nu
    sigma <- 1.0308 - 0.26758 * y
    z <- (w - mu)/sigma
    p <- 1 - pnorm(z)
    ret.list <- c(ret.list, list(statistic = W, parameters = n, 
        z.value = z, p.value = p, alternative = paste("True cdf does not equal the\n", 
            space(33), ret.list$distribution, " Distribution.", 
            sep = ""), method = "Shapiro-Francia GOF", data = x))
    names(ret.list$statistic) <- "W'"
    names(ret.list$parameters) <- "n"
    ret.list <- ret.list[c("distribution", "dist.abb", "distribution.parameters", 
        "n.param.est", "estimation.method", "statistic", "sample.size", 
        "parameters", "z.value", "p.value", "alternative", "method", 
        "data", "data.name", "bad.obs")]
    oldClass(ret.list) <- "gof"
    ret.list
}
