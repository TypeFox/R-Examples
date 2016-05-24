`print.crossval` <- function(x, digits = min(3, getOption("digits") - 3),
                             ...) {
    names(x$CVfit) <- paste("COCA", 1:x$n.axes, sep = "")
    cat("\nCross-validation for Predictive Co-Correspondence Analysis\n\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat(sprintf("\nCross-validatory %%fit of %s to %s:\n\n", x$nam.dat$namX,
                x$nam.dat$namY))
    print(round(x$CVfit, digits), print.gap = 2)
    invisible(x)
}

