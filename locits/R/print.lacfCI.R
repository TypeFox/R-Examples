print.lacfCI <-
function (x, ...) 
{
    cat("Class 'lacfCI' : Localized Autocovariance Object with Confidence Intervals:\n")
    cat("       ~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("Time localized to: ", x$nz, "\n")
    cat("\nsummary(.):\n----------\n")
    summary.lacfCI(x)
}
