print.estimate <-
function (x, conf.cov.sig.digits = .Options$digits, limits.sig.digits = .Options$digits, 
    ...) 
{
    coll.string <- paste("\n", space(33), sep = "")
    cat("\nResults of Distribution Parameter Estimation\n")
    cat("--------------------------------------------\n\n")
    cat("Assumed Distribution:", space(12), x$distribution, "\n\n", 
        sep = "")
    if (!is.null(x$par)) {
        cat("Estimated Parameter(s):", space(10), paste(paste(format(names(x$par), 
            justify = "left"), format(x$par, ..., nsmall = 0), 
            sep = " = "), collapse = paste("\n", space(33), sep = "")), 
            "\n\n", sep = "")
        if (!is.null(x$method)) 
            cat("Estimation Method:", space(15), x$method, "\n\n", 
                sep = "")
    }
    if (!is.null(x$quantiles)) {
        cat("Estimated Quantile(s):", space(11), paste(paste(format(names(x$quantiles), 
            justify = "left"), format(x$quantiles, ..., nsmall = 0), 
            sep = " = "), collapse = paste("\n", space(33), sep = "")), 
            "\n\n", sep = "")
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
    if (!is.null(x$sample.size)) {
        if (length(x$sample.size) > 1) {
            cat("Sample Sizes:", space(20), paste(paste(format(names(x$sample.size), 
                justify = "left"), format(x$sample.size, nsmall = 0, 
                ...), sep = " = "), collapse = coll.string), 
                "\n\n", sep = "")
        }
        else {
            cat("Sample Size:", space(21), x$sample.size, "\n\n", 
                sep = "")
        }
    }
    if (!is.null(x$bad.obs) && any(x$bad.obs > 0)) {
        if (length(x$bad.obs) > 1) 
            cat("Number NA/NaN/Inf's:", space(13), paste(paste(format(names(x$bad.obs), 
                justify = "left"), format(x$bad.obs, nsmall = 0, 
                ...), sep = " = "), collapse = coll.string), 
                "\n\n", sep = "")
        else cat("Number NA/NaN/Inf's:", space(13), x$bad.obs, 
            "\n\n", sep = "")
    }
    if (!is.null(x$interval)) {
        print.intervalEstimate(x$interval, conf.cov.sig.digits = conf.cov.sig.digits, 
            limits.sig.digits = limits.sig.digits, ...)
    }
    invisible(x)
}
