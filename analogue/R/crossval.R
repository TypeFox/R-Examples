## crossval() - crossvalidation methods for transfer functions
`crossval` <- function(obj, ...) {
    UseMethod("crossval")
}

`print.crossval` <- function(x, digits = min(getOption("digits"), 5), ...) {
    writeLines(strwrap("Model Cross-validation:", prefix = "\t"), sep = "\n\n")
    print(x$call)
    cat("\n")
    method <- x$CVparams$method
    writeLines(strwrap(paste(" Method:", method)))
    if(identical(method, "bootstrap"))
        writeLines(strwrap(paste(" No. Bootstraps:", x$CVparams$nboot)))
    if(identical(method, "kfold")) {
        writeLines(strwrap(paste("k:", x$CVparams$nfold)))
        writeLines(strwrap(paste("No. of folds:", x$CVparams$folds)))
    }
    cat("\n")
    print(x$performance, digits = digits, ...)
}
