varTest <-
function (x, alternative = "two.sided", conf.level = 0.95, sigma.squared = 1, 
    data.name = NULL) 
{
    if (is.null(data.name)) 
        data.name <- deparse(substitute(x))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if (!missing(conf.level)) 
        if ((length(conf.level) != 1) || !is.finite(conf.level) || 
            (conf.level <= 0) || (conf.level >= 1)) 
            stop("argument 'conf.level' must be a single number greater than zero and less than one.")
    alpha <- 1 - conf.level
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    num.df <- length(x) - 1
    if (length(unique(x)) < 2) {
        stop("All values in 'x' are the same.  Variance is zero.")
    }
    if (!missing(sigma.squared)) 
        if ((length(sigma.squared) != 1) || !is.finite(sigma.squared) || 
            sigma.squared <= 0) 
            stop(paste("argument 'sigma.squared' must be a", 
                "single positive numeric value."))
    var.x <- var(x)
    X <- (num.df * var.x)/sigma.squared
    p.less <- pchisq(X, df = num.df)
    p.greater <- 1 - pchisq(X, df = num.df)
    p.value <- switch(alternative, two.sided = 2 * min(p.less, 
        p.greater), less = p.less, greater = p.greater)
    statistic <- X
    names(statistic) <- "Chi-Squared"
    parameters <- num.df
    names(parameters) <- "df"
    null.value <- sigma.squared
    names(null.value) <- "variance"
    method <- "Chi-Squared Test on Variance"
    estimate <- var.x
    names(estimate) <- "variance"
    ci.type <- switch(alternative, two.sided = "two-sided", less = "upper", 
        greater = "lower")
    ci.interval <- enorm(x, ci = TRUE, ci.type = ci.type, conf.level = conf.level, 
        ci.param = "var")$interval$limits
    attr(ci.interval, "conf.level") <- conf.level
    ret.val <- c(list(statistic = statistic, parameters = parameters, 
        p.value = p.value, estimate = estimate, null.value = null.value, 
        alternative = alternative, method = method, data.name = data.name), 
        list(conf.int = ci.interval))
    oldClass(ret.val) <- "htest"
    return(ret.val)
}
