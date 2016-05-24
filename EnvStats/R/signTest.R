signTest <-
function (x, y = NULL, alternative = "two.sided", mu = 0, paired = FALSE, 
    conf.level = 0.95) 
{
    if (!is.null(y) && !paired) 
        stop(paste("Only the one-sample and paired sign test are available.", 
            "You must specify paired=T when", "supplying both 'x' and 'y'."))
    alternative <- match.arg(alternative, c("two.sided", "greater", 
        "less"))
    if (is.null(y)) {
        if (paired) 
            stop("argument 'y' missing for paired test.")
        if (!is.numeric(x)) 
            stop("'x' must be a numeric vector")
        data.name <- deparse(substitute(x))
        if (!is.numeric(mu) || length(mu) != 1 || !is.finite(mu)) 
            stop("argument 'mu' must be a single finite numeric value.")
        if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
            is.not.finite.warning(x)
            x <- x[x.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
        }
        x <- x[x != mu]
        n <- length(x)
        if (n < 2 || length(unique(x)) < 2) 
            stop(paste("'x' must contain at least 2 non-missing distinct values", 
                "not equal to the supplied value of 'mu'.", "This is not true for 'x' =", 
                data.name))
        d <- x
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
        }
        d <- d[d != mu]
        n <- length(d)
        if (n < 2 || all(d == d[1])) 
            stop(paste("There must be at least 2 non-missing distinct differences", 
                "not equal to the supplied value of 'mu'.", "This is not true for the paired values in\n", 
                "'x' =", data.name["x"], "and\n", "'y' =", data.name["y"]))
    }
    muhat <- median(d)
    names(muhat) <- ifelse(paired, "median of differences", "median")
    stat <- sum(d > mu)
    names(stat) <- ifelse(paired, "# Diffs > median of differences", 
        "# Obs > median")
    binom.list <- binom.test(x = stat, n = n, alternative = alternative)
    null.value <- mu
    names(null.value) <- ifelse(paired, "median of differences", 
        "median")
    ret.list <- list(statistic = stat, parameters = binom.list$parameters, 
        p.value = binom.list$p.value, estimate = muhat, null.value = null.value, 
        alternative = alternative, method = ifelse(paired, "Paired Sign test", 
            "Sign test"), data.name = data.name)
    ci.type <- switch(alternative, two.sided = "two-sided", greater = "lower", 
        less = "upper")
    ci.list <- eqnpar(x = d, p = 0.5, ci = TRUE, ci.type = ci.type, 
        approx.conf.level = conf.level)
    ci.list$interval$parameter <- ifelse(paired, "median of differences", 
        "median")
    ret.list <- c(ret.list, list(interval = ci.list$interval))
    oldClass(ret.list) <- "htest"
    ret.list
}
