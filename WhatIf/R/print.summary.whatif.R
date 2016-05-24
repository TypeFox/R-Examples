print.summary.whatif <- function(x, ...)  {

    x$sum.df$in.hull <- as.character(x$sum.df$in.hull)
    names(x$sum.df) <- c("Counterfactual", "In Hull", "Percent Nearby")
    
    #PRINT ON SCREEN
    cat("\nSummary of Counterfactual Inference Analysis\n")
    cat("\n")
    cat("Call:  ", deparse(x$call), "\n", sep="")
    cat("\n")
    cat("Total Number of Counterfactuals:  ", x$m, "\n", sep = "")
    cat("\n")
    cat("Number of Counterfactuals in Convex Hull:  ", x$m.inhull, "\n", sep = "")
    cat("\n")
    cat("Average Percent 'Nearby':  ", x$mean.near, "\n", sep = "")
    cat("\n")
    cat("Counterfactual in Convex Hull, True or False, and",
      " Percentage of Observed\n", sep = "")
    cat("Data Points 'Nearby' Counterfactual:\n")
    prmatrix(x$sum.df, rowlab = rep("", x$m), quote = FALSE)
    cat("\n")

    return(invisible(x))
}
