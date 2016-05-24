print.orest <- function (x, digits = max(3, getOption("digits") - 3), scientific=FALSE, ...) 
{
    cat("\nConstrained estimate:\n")
    if (length(x$b.restr))
        print.default(format(x$b.restr, digits = digits, scientific=scientific), print.gap = 2, 
            quote = FALSE, ...)
    else cat("No estimate available\n")
    cat("\n")
    invisible(x)
}
