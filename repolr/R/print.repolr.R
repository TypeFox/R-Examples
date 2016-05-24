print.repolr <-
function(x, digits = 4, ...){

    # output call and coefficients
    cat("\n")
    cat(x$title,"\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (length(x$coefficients)) {
        cat("Coefficients:", "\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2L, 
            quote = FALSE)
    }
    else cat("No Coefficients", "\n")
    cat("\n")
    invisible(x)
}
