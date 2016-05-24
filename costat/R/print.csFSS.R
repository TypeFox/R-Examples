print.csFSS <-
function (x, ...) 
{
    cat("Class 'csFSS' : Stationary Solutions Object from costat:\n")
    cat("       ~~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("\nsummary(.):\n----------\n")
    summary.csFSS(x)
}
