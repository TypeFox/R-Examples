chenTTest <-
function (x, y = NULL, alternative = "greater", mu = 0, paired = !is.null(y), 
    conf.level = 0.95, ci.method = "z") 
{
    if (!is.null(y) && !paired) 
        stop(paste("Only the one-sample and paired t-test are available", 
            "for Chen's test.  You must specify paired=T when", 
            "supplying both 'x' and 'y'."))
    alternative <- match.arg(alternative, c("greater", "less"))
    ci.method <- match.arg(ci.method, c("z", "t", "Avg. of z and t"))
    if (!missing(mu)) 
        if ((length(mu) != 1) || !is.finite(mu)) 
            stop("argument 'mu' must be a single finite numeric value.")
    skew.direction <- ifelse(alternative == "greater", "Positively", 
        "Negatively")
    if (is.null(y)) {
        if (paired) 
            stop("argument 'y' missing for paired test.")
        if (!is.vector(x, mode = "numeric") || is.factor(x)) 
            stop("'x' must be a numeric vector")
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
        muhat <- mean(x)
        sdhat <- sd(x)
        skewhat <- skewness(x, method = "fisher")
        stat.df.p.value.list <- chenTTest.sub(mu = mu, muhat = muhat, 
            sdhat = sdhat, skewhat = skewhat, n = n, alternative = alternative)
        method <- paste("One-sample t-Test\n", space(33), "Modified for\n", 
            space(33), skew.direction, "-Skewed Distributions\n", 
            space(33), "(Chen, 1995)", sep = "")
        ret.val <- c(stat.df.p.value.list, list(estimate = c(muhat, 
            sdhat, skewhat), null.value = mu, alternative = alternative, 
            method = method, sample.size = n, data.name = data.name, 
            bad.obs = bad.obs))
        names(ret.val$estimate) <- c("mean", "sd", "skew")
        names(ret.val$null.value) <- "mean"
    }
    else {
        if (!is.vector(x, mode = "numeric") || is.factor(x)) 
            stop("'x' must be a numeric vector")
        if (!is.vector(y, mode = "numeric") || is.factor(y)) 
            stop("'y' must be a numeric vector")
        if ((n <- length(x)) != length(y)) 
            stop("'x' and 'y' must have the same length when paired=TRUE.")
        data.name <- c(deparse(substitute(x)), deparse(substitute(y)))
        names(data.name) <- c("x", "y")
        d <- x - y
        if ((bad.obs <- sum(!(both.ok <- is.finite(d)))) > 0) {
            if (!all(is.finite(x))) 
                is.not.finite.warning(x)
            if (!all(is.finite(y))) 
                is.not.finite.warning(y)
            d <- d[both.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' or 'y' removed."))
            n <- length(d)
        }
        if (n < 2 || all(d == d[1])) 
            stop(paste("There must be at least 2 non-missing distinct differences. ", 
                "This is not true for the paired values in\n", 
                "'x' =", data.name["x"], "and\n", "'y' =", data.name["y"]))
        muhat <- mean(d)
        sdhat <- sd(d)
        skewhat <- skewness(d, method = "fisher")
        stat.df.p.value.list <- chenTTest.sub(mu = mu, muhat = muhat, 
            sdhat = sdhat, skewhat = skewhat, n = n, alternative = alternative)
        method <- paste("Paired-sample t-Test\n", space(33), 
            "Modified for\n", space(33), skew.direction, "-Skewed Distributions\n", 
            space(33), "(Chen, 1995)", sep = "")
        ret.val <- c(stat.df.p.value.list, list(estimate = c(muhat, 
            sdhat, skewhat), null.value = mu, alternative = alternative, 
            method = method, sample.size = n, data.name = data.name, 
            bad.obs = bad.obs))
        names(ret.val$estimate) <- c("mean of differences", "sd of differences", 
            "skew of differences")
        names(ret.val$null.value) <- "mean of differences"
    }
    ret.val <- ret.val[c("statistic", "parameters", "p.value", 
        "estimate", "null.value", "alternative", "method", "sample.size", 
        "data.name", "bad.obs")]
    ci.obj <- chenTTest.ci(muhat = muhat, sdhat = sdhat, skewhat = skewhat, 
        n = n, alternative = alternative, conf.level = conf.level, 
        p.value.type = ci.method, paired = paired)
    if (paired) 
        ci.obj$paramter <- "mean of differences"
    ret.val <- c(ret.val, list(interval = ci.obj))
    oldClass(ret.val) <- "htest"
    return(ret.val)
}
