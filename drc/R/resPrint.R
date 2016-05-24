"resPrint" <- function(resMat, headerText, interval, intervalLabel, display)
{
    if (display)
    {
        cat("\n")
        cat(paste(headerText, "\n", sep = ""))
        if (!identical(interval, "none"))
        {
            intervalText <- paste("(", intervalLabel, "-based confidence interval(s))\n", sep = "")
            cat(intervalText)
        }
        cat("\n")
        printCoefmat(resMat)
    }
#    invisible(resMat)
}