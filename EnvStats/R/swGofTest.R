swGofTest <-
function (x, distribution = c("norm", "lnorm", "lnormAlt", "lnorm3", 
    "zmnorm", "zmlnorm", "zmlnormAlt"), est.arg.list = NULL) 
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
    if (distribution == "lnorm3") {
        if (n < 5 || length(unique(x)) < 3) 
            stop(paste("'x' must contain at least 5 non-missing values,", 
                "and at least 3 distinct values."))
        if (n > 2000) 
            warning(paste("Too many observations.  This approximation only works", 
                "if the number of observations is between 5 and 2000"))
    }
    else {
        if (n < 4 || length(unique(x)) < 2) 
            stop(paste("'x' must contain at least 4 non-missing values,", 
                "and at least 2 distinct values."))
        if (n > 5000) 
            warning(paste("Too many observations.  This approximation only works", 
                "if the number of observations is between 4 and 5000"))
    }
    if (distribution == "lnorm3") 
        ret.list <- elnorm3(x, method = "zero.skew")
    else {
        est.fcn <- paste("e", distribution, sep = "")
        ret.list <- do.call(est.fcn, c(list(x = x), est.arg.list))
    }
    nrl <- names(ret.list)
    names(ret.list)[match("parameters", nrl)] <- "distribution.parameters"
    names(ret.list)[match("method", nrl)] <- "estimation.method"
    ret.list$data.name <- data.name
    ret.list$bad.obs <- bad.obs
    ret.list$dist.abb <- distribution
    new.x <- switch(distribution, norm = sort(x), lnorm = , lnormAlt = sort(log(x)), 
        lnorm3 = sort(log(x - ret.list$distribution.parameters["threshold"])), 
        zmnorm = sort(x[x != 0]), zmlnorm = , zmlnormAlt = sort(log(x[x > 
            0])))
    W <- swGofTestStatistic(new.x)
    if (n <= 11) {
        gam <- -2.273 + 0.459 * n
        w <- -log(gam - log(1 - W))
        mu <- 0.544 - 0.39978 * n + 0.025054 * n^2 - 0.0006714 * 
            n^3
        sigma <- exp(1.3822 - 0.77857 * n + 0.062767 * n^2 - 
            0.0020322 * n^3)
    }
    else {
        w <- log(1 - W)
        y <- log(n)
        mu <- -1.5861 - 0.31082 * y - 0.083751 * y^2 + 0.0038915 * 
            y^3
        sigma <- exp(-0.4803 - 0.082676 * y + 0.0030302 * y^2)
    }
    z <- (w - mu)/sigma
    if (distribution == "lnorm3") {
        u <- log(n)
        sdlog <- ret.list$distribution.parameters["sdlog"]
        v <- u * (sdlog - sdlog^2)
        if (n <= 11) {
            mu.z <- -3.8267 + 2.8242 * u - 0.63673 * u^2 - 0.020815 * 
                v
            sigma.z <- -4.9914 + 8.6724 * u - 4.27905 * u^2 + 
                0.7035 * u^3 - 0.013431 * v
        }
        else {
            mu.z <- -3.7796 + 2.4038 * u - 0.6675 * u^2 + 0.082863 * 
                u^3 - 0.0037935 * u^4 - 0.027027 * v - 0.0019887 * 
                v * u
            sigma.z <- 2.1924 - 1.0957 * u + 0.33737 * u^2 - 
                0.043201 * u^3 + 0.0019974 * u^4 - 0.0053312 * 
                v * u
        }
        z <- (z - mu.z)/sigma.z
    }
    p <- 1 - pnorm(z)
    ret.list <- c(ret.list, list(statistic = W, parameters = n, 
        z.value = z, p.value = p, alternative = paste("True cdf does not equal the\n", 
            space(33), ret.list$distribution, " Distribution.", 
            sep = ""), method = "Shapiro-Wilk GOF", data = x))
    names(ret.list$statistic) <- "W"
    names(ret.list$parameters) <- "n"
    ret.list <- ret.list[c("distribution", "dist.abb", "distribution.parameters", 
        "n.param.est", "estimation.method", "statistic", "sample.size", 
        "parameters", "z.value", "p.value", "alternative", "method", 
        "data", "data.name", "bad.obs")]
    oldClass(ret.list) <- "gof"
    ret.list
}
