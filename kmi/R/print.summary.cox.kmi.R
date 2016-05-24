print.summary.cox.kmi <- function(x, digits = max(getOption("digits") - 3, 3),
                                  signif.stars = getOption("show.signif.stars"),
                                  print.ind = FALSE, ...) {
    if (!inherits(x, "summary.cox.kmi")) {
        stop("'x' must be of class 'kmi'")
    }
    cat("Call:\n")
    dput(x$call)
    cat("\n")
    cat("\n")
    savedig <- options(digits = digits)
    on.exit(options(savedig))
    cat("*****************\n")
    cat("Pooled estimates:\n")
    cat("*****************\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars=signif.stars, ...)
    cat("\n")
    print(x$conf.int)
    cat("\n")
    if (print.ind) {
        cat("*********************\n")
        cat("Individual estimates:\n")
        cat("*********************\n\n")
        for (i in seq_along(x$individual.fit)) {
            cat(paste("*** Imputation", i, "***", sep = " ")); cat("\n")
            print(x$individual.fit[[i]], digits = max(getOption("digits") - 3, 3), ...)
            cat("\n")
        }
    }
    invisible()
}
