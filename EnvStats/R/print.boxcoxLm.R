print.boxcoxLm <-
function (x, ...) 
{
    cat("\nResults of Box-Cox Transformation\n")
    cat("---------------------------------\n\n")
    cat("Objective Name:", space(18), x$objective.name, "\n\n", 
        sep = "")
    data.name <- x$data.name
    cat("Linear Model:", space(20), data.name, "\n\n", sep = "")
    cat("Sample Size:", space(21), x$sample.size, "\n\n", sep = "")
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
