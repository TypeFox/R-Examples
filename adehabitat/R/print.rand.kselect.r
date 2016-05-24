"print.rand.kselect" <- function(x, ...)
{
    cat("****** Randomization tests of the k-select analysis ******\n\n")
    cat("Test of the first eigenvalue:\n")
    cat("Observed value:", x$global[1], "\n")
    cat("P-value :", x$global[2], "\n")
    cat("\n")
    cat("Test of the marginality of each individual\n(to be compared with bonferroni alpha level:",
        x$alpha/nrow(x$marg),"):\n\n")
    print(x$marg, ...)
    cat("\nSign of the mean for each animal and each variable: \n")
    cat("(when significant, the sign is tripled) \n\n")
    print(x$per.ind$signification, quote=FALSE)
    cat("\n\nOther elements of the list $per.ind:")
    cat("\n  $obsval: mean of variables for each animal")
    cat("\n  $pvalue: P-value of the means in $obsval\n\n")
}

