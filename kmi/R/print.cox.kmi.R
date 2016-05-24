print.cox.kmi <- function(x, print.ind = FALSE, ...) {
    if (!inherits(x, "cox.kmi")) {
        stop("'x' must an object of class 'kmi'")
    }
    if ("coxph.penal" %in% class(x$cox.kmi.fit[[1]])) {
        stop("Some things to think about before using frailties here")
    }
    cat("Call:\n")
    dput(x$call); cat("\n")
    coef <- x$coefficients
    se <- sqrt(diag(x$variance))
    tmp <- cbind(coef, exp(coef), se, x$coefficients / sqrt(diag(x$variance)),
                 2 * pt(abs(x$coefficients / sqrt(diag(x$variance))), df = x$df, lower.tail = FALSE))
    dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", "se(coef)", "t", "p"))
    cat("*****************\n")
    cat("Pooled estimates:\n")
    cat("*****************\n")
    print(tmp); cat("\n")
    if (print.ind) {
        cat("*********************\n")
        cat("Individual estimates:\n")
        cat("*********************\n\n")
        for (i in seq_along(x$individual.fit)) {
            cat(paste("*** Imputation", i, "***", sep = " ")); cat("\n")
            print(x$individual.fit[[i]], ...)
            cat("\n")
        }
    }
    invisible()
}
