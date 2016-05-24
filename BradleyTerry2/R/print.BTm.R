print.BTm <- function (x, ...)
{
    cat("Bradley Terry model fit by ")
    cat(x$method, "\n")
    NextMethod()
}
