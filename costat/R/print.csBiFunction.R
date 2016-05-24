print.csBiFunction <-
function (x, ...) 
{
    cat("Class 'csBiFunction' : Contains two sampled functions:\n")
    cat("       ~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("\nsummary(.):\n----------\n")
    summary.csBiFunction(x)
}
