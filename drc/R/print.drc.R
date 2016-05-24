"print.drc" <- function(x, ..., digits = max(3, getOption("digits") - 3)) 
{
    object <- x

    classList <- class(object)
    cat(paste("\n", "A 'drc' model.", "\n", sep=""))
    
    ## Borrowing from print.lm
    cat("\nCall:\n", deparse(object$"call"), "\n\n", sep = "")
    if (length(coef(object))>0) 
    {
        cat("Coefficients:\n")
        print.default(format(coef(object), digits = digits), print.gap = 2, quote = FALSE)
    } else {
        cat("No coefficients\n")
    }
    cat("\n")
    
    invisible(object)
}
