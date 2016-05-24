print.estimateCensored <-
function (x, show.cen.levels = TRUE, pct.censored.digits = .Options$digits, 
    conf.cov.sig.digits = .Options$digits, limits.sig.digits = .Options$digits, 
    ...) 
{
    coll.string <- paste("\n", space(33), sep = "")
    cat("\nResults of Distribution Parameter Estimation\n")
    cat("Based on Type I Censored Data\n")
    cat("--------------------------------------------\n\n")
    cat("Assumed Distribution:", space(12), x$distribution, "\n\n", 
        sep = "")
    cat("Censoring Side:", space(18), x$censoring.side, "\n\n", 
        sep = "")
    if (show.cen.levels) {
        cat("Censoring Level(s):", space(12), format(x$censoring.levels, 
            nsmall = 0, justify = "left", ...), "\n\n")
    }
    if (!is.null(x$parameters)) {
        cat("Estimated Parameter(s):", space(10), paste(paste(format(names(x$parameters), 
            justify = "left"), format(x$parameters, ..., nsmall = 0), 
            sep = " = "), collapse = coll.string), "\n\n", sep = "")
        if (!is.null(x$method)) 
            cat("Estimation Method:", space(15), x$method, "\n\n", 
                sep = "")
        if (!is.null(x$prob.method)) 
            cat("Plotting Position Method:", space(8), x$prob.method, 
                "\n\n", sep = "")
        if (!is.null(x$plot.pos.con)) 
            cat("Plotting Position Constant:", space(6), x$plot.pos.con, 
                "\n\n", sep = "")
    }
    if (!is.null(x$quantiles)) {
        cat("Estimated Quantile(s):", space(11), paste(paste(format(names(x$quantiles), 
            justify = "left"), format(x$quantiles, ..., nsmall = 0), 
            sep = " = "), collapse = coll.string), "\n\n", sep = "")
        if (!is.null(x$quantile.method)) 
            cat("Quantile Estimation Method:", space(6), x$quantile.method, 
                "\n\n", sep = "")
    }
    if (is.null(names(x$data.name))) 
        cat("Data:", space(28), x$data.name, "\n\n", sep = "")
    else cat("Data:", space(28), paste(paste(format(names(x$data.name), 
        justify = "left"), format(x$data.name, ...), sep = " = "), 
        collapse = coll.string), "\n\n", sep = "")
    if (!is.null(x$subset.expression)) 
        cat("Subset With:", space(21), x$subset.expression, "\n\n", 
            sep = "")
    if (!is.null(x$parent.of.data)) 
        cat("Data Source:", space(21), x$parent.of.data, "\n\n", 
            sep = "")
    cat("Censoring Variable:", space(14), x$censoring.name, "\n\n", 
        sep = "")
    if (!is.null(x$bad.obs) && x$bad.obs > 0) 
        cat("Number NA/NaN/Inf's Removed:", space(5), x$bad.obs, 
            "\n\n", sep = "")
    cat("Sample Size:", space(21), x$sample.size, "\n\n", sep = "")
    cat("Percent Censored:", space(16), round(x$percent.censored, 
        pct.censored.digits), "%", "\n\n", sep = "")
    if (!is.null(x$interval)) {
        print.intervalEstimateCensored(x$interval, conf.cov.sig.digits = conf.cov.sig.digits, 
            limits.sig.digits = limits.sig.digits, ...)
    }
    invisible(x)
}
