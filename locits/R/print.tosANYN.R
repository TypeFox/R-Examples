print.tosANYN <-
function (x, ...) 
{
    cat("Class 'tosANYN' : Stationarity Object for Arbitrary Length Data :\n")
    cat("       ~~~~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("\nsummary(.):\n----------\n")
    summary.tosANYN(x)
}
