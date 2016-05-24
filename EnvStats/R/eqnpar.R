eqnpar <-
function (x, p = 0.5, ci = FALSE, lcl.rank = NULL, ucl.rank = NULL, 
    lb = -Inf, ub = Inf, ci.type = "two-sided", ci.method = "exact", 
    approx.conf.level = 0.95, digits = 0) 
{
    if (!is.vector(x, mode = "numeric") & !is.factor(x)) 
        stop("'x' must be a numeric vector or a factor.")
    if (is.factor(x)) {
        x <- as.numeric(x)
        warning("x is a factor; numeric values of levels used.\n")
    }
    if (!is.vector(p, mode = "numeric")) 
        stop("'p' must be a numeric vector.")
    if (any(!is.finite(p))) 
        stop("NA/NaN/Inf values not allowed for 'p'")
    if (any(p < 0) || any(p > 1)) 
        stop("All values of 'p' must be between 0 and 1.")
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
    q <- quantile(x, p)
    if (length(p) == 1 && p == 0.5) 
        names(q) <- "Median"
    else {
        pct <- round(100 * p, digits)
        names(q) <- paste(pct, number.suffix(pct), " %ile", sep = "")
    }
    ret.list <- list(distribution = "None", sample.size = n, 
        method = "Nonparametric", quantiles = q, quantile.method = "Nonparametric", 
        data.name = data.name, bad.obs = bad.obs)
    if (ci) {
        if (any(length.list(lb, ub) != 1) || lb >= ub) 
            stop(paste("'lb' and 'ub' must be scalars,", "and 'lb' must be strictly less than 'ub'"))
        if (length(p) > 1) 
            stop("When 'ci' = TRUE, 'p' must be a scalar")
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        ci.method <- match.arg(ci.method, c("exact", "normal.approx"))
        if (!is.numeric(approx.conf.level) || length(approx.conf.level) != 
            1 || approx.conf.level <= 0 || approx.conf.level >= 
            1) 
            stop("'approx.conf.level' must be a scalar between 0 and 1")
        ci.obj <- ci.qnpar(x = x, p = p, lcl.rank = lcl.rank, 
            ucl.rank = ucl.rank, lb = lb, ub = ub, ci.type = ci.type, 
            ci.method = ci.method, approx.alpha = 1 - approx.conf.level, 
            digits = digits)
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    oldClass(ret.list) <- "estimate"
    ret.list
}
