"print.grouped" <-
function(x, ...){
    if(!inherits(x, "grouped"))
        stop("Use only with 'grouped' objects.\n")
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if(length(coefs <- x$coef)){
        coefs <- coefs[-length(coefs)]
        cat("Coefficients:\n")
        print.default(format(coefs, digits = 3), print.gap = 2, quote = FALSE)
    } else cat("No coefficients\n")
    cat("\n\n")
    invisible(x)
}

