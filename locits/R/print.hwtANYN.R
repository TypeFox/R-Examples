print.hwtANYN <-
function (x, ...) 
{
    cat("Class 'hwtANYN' : Haar Wavelet for Arbitrary Length Data object:\n")
    cat("       ~~~~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("\nsummary(.):\n----------\n")
    summary.hwtANYN(x)
}
