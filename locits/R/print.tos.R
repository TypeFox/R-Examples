print.tos <-
function (x, ...) 
{
    cat("Class 'tos' : Stationarity Object :\n")
    cat("       ~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("\nsummary(.):\n----------\n")
    summary.tos(x)
}
