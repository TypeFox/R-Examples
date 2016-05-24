##' @S3method print polywog
print.polywog <- function(x, ...)
{
    ## Print the function call used to fit the model (using same code as in
    ## 'print.lm')
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")

    ## Extract lists of expanded and non-expanded terms
    expanded <- attr(x$varNames, "expanded")
    inTerms <- paste(x$varNames[expanded], collapse = ", ")
    cat("Variables included in polynomial expansion of degree ", x$degree, ":\n",
        sep = "")
    writeLines(strwrap(inTerms, prefix = "  "))
    if (any(!expanded)) {
        outTerms <- paste(x$varNames[!expanded], collapse = ", ")
        cat("\nVariables included linearly:\n")
        writeLines(strwrap(outTerms, prefix = "  "))
    }

    ## Print regularization method
    cat("\nRegularization method:",
        switch(x$method, alasso = "Adaptive LASSO", scad = "SCAD"))

    cat("\n\n")
    invisible(x)
}
