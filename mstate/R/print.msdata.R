print.msdata <- function(x,trans=FALSE,...)
{
    if (!inherits(x, "msdata"))
        stop("'x' must be an 'msdata' object")
    cat("An object of class 'msdata'\n\nData:\n")
    print.data.frame(x)
    if (trans) {
        trans <- attr(x, "trans")
        cat("\nTransition matrix:\n")
        print(trans)
    }
    return(invisible())
}
