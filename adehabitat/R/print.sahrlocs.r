"print.sahrlocs" <- function(x, ...)
{
    ## Verifications
    if (!inherits(x, "sahrlocs")) stop("object should be of type \"sahrlocs\"")

    ## The output
    cat("************** Object of type sahrlocs **************\n\n")
    nr<-attr(x, "nrow")
    nc<-attr(x, "ncol")
    cat("The area of interest is a ", nr, "*", nc, " raster matrix\n")
    nc<-ncol(x$locs)
    cat(nc, " animals are available :\n")
    print(names(as.data.frame(unclass(x$hr))), ...)
    cat("\n\n the following variables are available for the study area:\n")
    print(names(as.data.frame(unclass(x$sa))), ...)

    if (!is.null(x$descan)) {
        cat("\nthe following variables are available for each monitored animal:\n")
        print(names(x$descan), ...)
    } else {
        cat("\nno variables have been measured on the animals\n")
    }

}

