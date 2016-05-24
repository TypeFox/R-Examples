print.boxcoxCensored <-
function (x, ...) 
{
    coll.string <- paste("\n", space(33), sep = "")
    cat("\nResults of Box-Cox Transformation\n")
    cat("Based on Type I Censored Data\n")
    cat("---------------------------------\n\n")
    cat("Objective Name:", space(18), x$objective.name, "\n\n", 
        sep = "")
    if (is.null(names(x$data.name))) 
        cat("Data:", space(28), x$data.name, "\n\n", sep = "")
    else cat("Data:", space(28), paste(paste(format(names(x$data.name), 
        justify = "left"), format(x$data.name, ...), sep = " = "), 
        collapse = coll.string), "\n\n", sep = "")
    if (!is.null(x$subset.expression)) 
        cat("Subset With:", space(21), x$subset.expression, "\n\n", 
            sep = "")
    cat("Censoring Variable:", space(14), x$censoring.name, "\n\n", 
        sep = "")
    cat("Censoring Side:", space(18), x$censoring.side, "\n\n", 
        sep = "")
    cat("Censoring Level(s):", space(12), format(x$censoring.levels, 
        nsmall = 0, justify = "left", ...), "\n\n")
    if (!is.null(x$parent.of.data)) 
        cat("Data Source:", space(21), x$parent.of.data, "\n\n", 
            sep = "")
    if (!is.null(x$bad.obs) && any(x$bad.obs > 0)) 
        cat("Number NA/NaN/Inf's Removed:", space(5), x$bad.obs, 
            "\n\n", sep = "")
    cat("Sample Size:", space(21), x$sample.size, "\n\n", sep = "")
    cat("Percent Censored:", space(16), round(x$percent.censored, 
        1), "%", "\n\n", sep = "")
    if (x$optimize) {
        cat("Bounds for Optimization:", space(9), paste(paste(format(names(x$optimize.bounds), 
            justify = "left"), format(x$optimize.bounds, nsmall = 0, 
            ...), sep = " = "), collapse = paste("\n", space(33), 
            sep = "")), "\n\n", sep = "")
        cat("Optimal Value:", space(19), paste("lambda =", format(x$lambda, 
            ...)), "\n\n", sep = "")
        cat("Value of Objective:", space(14), paste(x$objective.name, 
            "=", format(x$objective, ...)), "\n\n", sep = "")
    }
    else {
        dum.mat <- cbind(x$lambda, x$objective)
        dimnames(dum.mat) <- list(rep("", nrow(dum.mat)), c("lambda", 
            x$objective.name))
        print(dum.mat, ...)
    }
    invisible(x)
}
