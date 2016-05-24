print.lacf <-
function (x, ...) 
{
    cat("Class 'lacf' : Localized Autocovariance/correlation Object:\n")
    cat("       ~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("\nsummary(.):\n----------\n")
    summary.lacf(x)
}
