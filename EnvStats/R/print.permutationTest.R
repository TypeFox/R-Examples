print.permutationTest <-
function (x, ...) 
{
    coll.string <- paste("\n", space(33), sep = "")
    cat("\nResults of Hypothesis Test\n")
    cat("--------------------------\n\n")
    alt.string <- x$alternative
    if (!is.null(nv <- x$null.value)) {
        nnv <- names(nv)
        if ((lnv <- length(nv)) == 1) {
            cat("Null Hypothesis:", space(17), paste(paste(format(nnv, 
                justify = "left"), format(nv, nsmall = 0, ...), 
                sep = " = "), collapse = coll.string), "\n\n", 
                sep = "")
            alt.string <- x$alternative
            if (!is.na(match(alt.string, c("two.sided", "less", 
                "greater")))) {
                alt.string <- switch(alt.string, two.sided = "not equal to", 
                  less = "less than", greater = "greater than")
                alt.string <- paste("True", nnv, "is", alt.string, 
                  format(nv, nsmall = 0, ...))
            }
        }
        else {
            cat("Null Hypothesis:", space(17), paste("All", lnv, 
                "values of", nnv[1], "=", format(nv[1], ...)), 
                "\n\n", sep = "")
        }
    }
    cat("Alternative Hypothesis:", space(10), alt.string, "\n\n", 
        sep = "")
    cat("Test Name:", space(23), x$method, "\n\n", sep = "")
    if (!is.null(x$estimate)) {
        cat("Estimated Parameter(s):", space(10), paste(paste(format(names(x$estimate), 
            justify = "left"), format(x$estimate, nsmall = 0, 
            ...), sep = " = "), collapse = coll.string), "\n\n", 
            sep = "")
        if (!is.null(x$estimation.method)) 
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
    if (!is.null(x$statistic)) {
        string <- ifelse(length(x$statistic) == 1, paste("Test Statistic:", 
            space(18), sep = ""), paste("Test Statistics:", space(17), 
            sep = ""))
        cat(string, paste(paste(format(names(x$statistic), justify = "left"), 
            format(x$statistic, nsmall = 0, ...), sep = " = "), 
            collapse = coll.string), "\n\n", sep = "")
    }
    if (!is.null(x$parameters)) {
        string <- ifelse(length(x$parameters) > 1, paste("Test Statistic Parameters:", 
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
    if (!is.null(x$conf.int)) {
        ci.pct <- format(100 * attr(x$conf.int, "conf.level"))
        cat(ci.pct, "% Confidence Interval:", space(33 - nchar(ci.pct) - 
            22), paste(paste(c("LCL", "UCL"), format(x$conf.int, 
            nsmall = 0, ...), sep = " = "), collapse = coll.string), 
            "\n\n", sep = "")
    }
    if (!is.null(x$interval)) {
        print.intervalEstimate(x$interval)
    }
    invisible(x)
}
