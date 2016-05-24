summary.wt.filter <- function(object, ...)
{
    cat("Filter Class: ")
    cat(object@wt.class)
    cat("\n")
    cat("Name: ")
    cat(toupper(object@wt.name))
    cat("\n")
    cat("Length: ")
    cat(object@L)
    cat("\n")

    invisible(object)
}

