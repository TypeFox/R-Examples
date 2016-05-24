skewGofTest <-
function (x, distribution = c("norm", "lnorm", "lnormAlt", "zmnorm", 
    "zmlnorm", "zmlnormAlt"), est.arg.list = NULL, alternative = c("two.sided", 
    "less", "greater")) 
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
    if (n < 8) 
        stop("'x' must have at least 8 non-missing values")
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
    sqrt.beta1 <- skewness(new.x, method = "moment")
    alternative <- match.arg(alternative)
    z <- sqrt.beta1/sqrt((6 * (n - 2))/((n + 1) * (n + 3)))
    if (n < 150) {
        beta2 <- (3 * (n^2 + 27 * n - 70) * (n + 1) * (n + 3))/((n - 
            2) * (n + 5) * (n + 7) * (n + 9))
        W <- sqrt(-1 + sqrt(2 * beta2 - 2))
        delta <- 1/sqrt(log(W))
        alpha <- sqrt(2/(W^2 - 1))
        z <- z/alpha
        z <- delta * log(z + sqrt(z^2 + 1))
    }
    p <- switch(alternative, two.sided = 2 * pnorm(-abs(z)), 
        less = pnorm(z), greater = 1 - pnorm(z))
    alt.string <- switch(distribution, norm = "True skew ", lnorm = , 
        lnormAlt = paste("True skew of log-transformed\n", space(33), 
            "values ", sep = ""), zmnorm = "True skew of non-zero values ", 
        zmlnorm = , zmlnormAlt = paste("True skew of log-transformed\n", 
            space(33), "positive values ", sep = ""))
    ret.list <- c(ret.list, list(statistic = sqrt.beta1, parameters = n, 
        z.value = z, p.value = p, alternative = paste(alt.string, 
            ifelse(alternative == "two.sided", "does not equal ", 
                paste("is", alternative, "than ")), "0.", sep = ""), 
        method = "Zero-Skew GOF", data = x))
    names(ret.list$statistic) <- switch(distribution, norm = "Skew", 
        lnorm = , lnormAlt = "Skew[Log(x)]", zmnorm = "Skew(x|x!=0)", 
        zmlnorm = , zmlnormAlt = "Skew[ Log(x|x>0) ]")
    names(ret.list$parameters) <- "n"
    ret.list <- ret.list[c("distribution", "dist.abb", "distribution.parameters", 
        "n.param.est", "estimation.method", "statistic", "sample.size", 
        "parameters", "z.value", "p.value", "alternative", "method", 
        "data", "data.name", "bad.obs")]
    oldClass(ret.list) <- "gof"
    ret.list
}
