print.lacv <-
function (x, ...) 
{
    cat("Class 'lacv' : Localized Autocovariance/correlation Object:\n")
    cat("       ~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("\nsummary(.):\n----------\n")
    summary.lacv(x)
}
