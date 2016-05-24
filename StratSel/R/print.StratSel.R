print.StratSel <-
function(x, ...)
{
    cat("Strategic Selection Model: StratSel \n")
    cat(" \n")
    cat("Number of Obs\n")
    print(x$df)
    cat("\n Coefficients:\n")
    print.summary.StratSel(x$coefficients)
}
