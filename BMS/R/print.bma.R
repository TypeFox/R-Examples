print.bma <-
function (x, ...) 
{
    if (!is.bma(x)) {
        return(print(x))
    }
    print(estimates.bma(x), include.constant = TRUE, ...)
    cat("\n")
    print(info.bma(x), ...)
    cat("\n")
}
