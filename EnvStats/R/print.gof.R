print.gof <-
function (x, ...) 
{
    coll.string <- paste("\n", space(33), sep = "")
    cat("\nResults of Goodness-of-Fit Test\n")
    cat("-------------------------------\n\n")
    cat("Test Method:", space(21), x$method, "\n\n", sep = "")
    cat("Hypothesized Distribution:", space(7), x$distribution, 
        "\n\n", sep = "")
    if (!is.null(x$n.param.est) && x$n.param.est > 0) {
        cat("Estimated Parameter(s):", space(10), paste(paste(format(names(x$distribution.parameters), 
            justify = "left"), format(x$distribution.parameters, 
            nsmall = 0, ...), sep = " = "), collapse = coll.string), 
            "\n\n", sep = "")
        cat("Estimation Method:", space(15), x$estimation.method, 
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
    if (!is.null(x$bad.obs) && any(x$bad.obs > 0)) {
        if (length(x$bad.obs) == 1) 
            cat("Number NA/NaN/Inf's Removed:", space(5), x$bad.obs, 
                "\n\n", sep = "")
        else {
            cat("Number NA/NaN/Inf's Removed:", space(5), paste(paste(format(names(x$bad.obs), 
                justify = "left"), format(x$bad.obs, nsmall = 0, 
                ...), sep = " = "), collapse = coll.string), 
                "\n\n", sep = "")
        }
    }
    if (length(x$sample.size) > 1) {
        cat("Sample Sizes:", space(20), paste(paste(format(names(x$sample.size), 
            justify = "left"), format(x$sample.size, nsmall = 0, 
            ...), sep = " = "), collapse = coll.string), "\n\n", 
            sep = "")
    }
    else {
        cat("Sample Size:", space(21), x$sample.size, "\n\n", 
            sep = "")
    }
    string <- ifelse(length(x$statistic) == 1, paste("Test Statistic:", 
        space(18), sep = ""), paste("Test Statistics:", space(17), 
        sep = ""))
    cat(string, paste(paste(format(names(x$statistic), justify = "left"), 
        format(x$statistic, nsmall = 0, ...), sep = " = "), collapse = coll.string), 
        "\n\n", sep = "")
    n.params <- length(x$parameters)
    if (n.params > 0) {
        string <- ifelse(n.params > 1, paste("Test Statistic Parameters:", 
            space(7), sep = ""), paste("Test Statistic Parameter:", 
            space(8), sep = ""))
        cat(string, paste(paste(format(names(x$parameters), justify = "left"), 
            format(x$parameters, nsmall = 0, ...), sep = " = "), 
            collapse = coll.string), "\n\n", sep = "")
    }
    if (length(x$p.value) == 1) 
        cat("P-value:", space(25), format(x$p.value, ...), "\n\n", 
            sep = "")
    else {
        if (!is.null(names(x$p.value))) 
            cat("P-values:", space(24), paste(paste(format(names(x$p.value), 
                justify = "left"), format(x$p.value, ...), sep = " = "), 
                collapse = coll.string), "\n\n", sep = "")
        else cat("P-values:", space(24), paste(format(x$p.value, 
            ...), collapse = coll.string), "\n\n", sep = "")
    }
    cat("Alternative Hypothesis:", space(10), x$alternative, 
        "\n", sep = "")
    invisible(x)
}
