print.BootStep <-
function (x, ...) {
    cat("\nSummary of Bootstrapping the 'stepAIC()' procedure for\n")
    cat("\nCall:\n", paste(deparse(x$OrigModel$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat("Bootstrap samples:", length(x$BootStepAIC), "\n")
    cat("Direction:", x$direction, "\n")
    cat("Penalty:", round(x$k, 2), "* df\n")
    cat("\nCovariates selected\n")
    print(round(x$Covariates, 2))
    cat("\nCoefficients Sign\n")
    print(round(x$Sign, 2))
    cat("\nStat Significance\n")
    print(round(x$Significance, 2))    
    cat("\n\nThe stepAIC() for the original data-set gave\n")
    if (inherits(x$OrigModel, "polr"))
        cat("\n")
    print(x$OrigStepAI)
    cat("\n")
    print(x$OrigStepAI$anova)
    cat("\n")
    invisible(x)    
}

