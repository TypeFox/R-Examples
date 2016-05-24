print.htestCensored <-
function (x, show.cen.levels = TRUE, ...) 
{
    coll.string <- paste("\n", space(33), sep = "")
    cat("\nResults of Hypothesis Test\n")
    cat("Based on Censored Data\n")
    cat("--------------------------\n\n")
    alt.string <- x$alternative
    if (!is.null(nv <- x$null.value)) {
        nnv <- names(nv)
        if ((lnv <- length(nv)) == 1) {
            cat("Null Hypothesis:", space(17), paste(paste(format(nnv, 
                justify = "left"), format(nv, ...), sep = " = "), 
                collapse = coll.string), "\n\n", sep = "")
            alt.string <- x$alternative
            if (!is.na(match(alt.string, c("two.sided", "less", 
                "greater")))) {
                alt.string <- switch(alt.string, two.sided = "not equal to", 
                  less = "less than", greater = "greater than")
                alt.string <- paste("True", nnv, "is", alt.string, 
                  nv)
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
    cat("Censoring Side:", space(18), x$censoring.side, "\n\n", 
        sep = "")
    if (show.cen.levels) {
        if ((n <- length(x$censoring.levels)) == 1) 
            cat("Censoring Level(s):", space(12), format(x$censoring.levels, 
                nsmall = 0, ...), "\n\n")
        else {
            ncl <- names(x$censoring.levels)
            cat("Censoring Level(s):", space(12), ncl[1], "=", 
                format(x$censoring.levels[[1]], nsmall = 0, ...), 
                "\n")
            for (i in 2:n) cat(space(32), ncl[i], "=", format(x$censoring.levels[[i]], 
                nsmall = 0, ...), "\n")
            cat("\n")
        }
    }
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
    if (is.null(names(x$censoring.name))) 
        cat("Censoring Variable:", space(14), x$censoring.name, 
            "\n\n", sep = "")
    else cat("Censoring Variable:", space(14), paste(paste(format(names(x$censoring.name), 
        justify = "left"), format(x$censoring.name, ...), sep = " = "), 
        collapse = coll.string), "\n\n", sep = "")
    if (!is.null(x$bad.obs) && any(x$bad.obs > 0)) {
        if (length(x$bad.obs) > 1) 
            cat("Number NA/NaN/Inf's Removed:", space(5), paste(paste(format(names(x$bad.obs), 
                justify = "left"), format(x$bad.obs, nsmall = 0, 
                ...), sep = " = "), collapse = coll.string), 
                "\n\n", sep = "")
        else cat("Number NA/NaN/Inf's Removed:", space(5), x$bad.obs, 
            "\n\n", sep = "")
    }
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
    if (is.null(names(x$percent.censored))) 
        cat("Percent Censored:", space(16), round(x$percent.censored, 
            1), "%", "\n\n", sep = "")
    else cat("Percent Censored:", space(16), paste(paste(format(names(x$percent.censored), 
        justify = "left"), paste(format(round(x$percent.censored, 
        1), nsmall = 0, ...), "%", sep = ""), sep = " = "), collapse = coll.string), 
        "\n\n", sep = "")
    string <- ifelse(length(x$statistic) == 1, paste("Test Statistic:", 
        space(18), sep = ""), paste("Test Statistics:", space(17), 
        sep = ""))
    cat(string, paste(paste(format(names(x$statistic), justify = "left"), 
        format(x$statistic, nsmall = 0, ...), sep = " = "), collapse = coll.string), 
        "\n\n", sep = "")
    if (!is.null(x$parameter)) {
        string <- ifelse(length(x$parameter) > 1, paste("Test Statistic Parameters:", 
            space(7), sep = ""), paste("Test Statistic Parameter:", 
            space(8), sep = ""))
        cat(string, paste(paste(format(names(x$parameter), justify = "left"), 
            format(x$parameter, nsmall = 0, ...), sep = " = "), 
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
        ci.string <- format(100 * attr(x$conf.int, "conf.level"))
        ncci <- nchar(ci.string)
        cat(ci.string, "% Confidence Interval:", space(33 - ncci - 
            22), paste(paste(c("LCL", "UCL"), format(x$conf.int, 
            nsmall = 0, ...), sep = " = "), collapse = coll.string), 
            "\n\n", sep = "")
    }
    invisible(x)
}
