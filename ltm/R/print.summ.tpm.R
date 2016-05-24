print.summ.tpm <-
function (x, digits = max(3, getOption("digits") - 3), ...){
    if (!inherits(x, "summ.tpm"))
        stop("Use only with 'summ.tpm' objects.\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat("Model Summary:\n")
    model.sum <- data.frame(log.Lik = x$logLik, AIC = x$AIC, BIC = x$BIC, row.names = "")
    print(model.sum)
    cat("\nCoefficients:\n")
    if (is.null(x$attr)) {
        print(round(x$coef, digits = digits))
    } else {
        dat <- data.frame(x$coef, check.names = FALSE)
        dat[] <- lapply(dat, round, digits = digits)
        nr <- nrow(dat)
        dat$" " <- if (x$type == "rasch") c(rep(x$attr, length.out = nr - 1), "") else rep(x$attr, length.out = nr)
        print(dat)
        if (x$ancr) cat("\n'*' denotes an anchoring item\n\n") else cat("\n'*' denotes a linking item\n\n")
    }
    cat("\nIntegration:\n")
    cat("method: Gauss-Hermite\n")
    cat("quadrature points:", x$control$GHk, "\n")
    cat("\nOptimization:\n")
    cat("Optimizer:", if (x$control$optimizer == "optim") paste("optim (", x$control$method, ")\n", sep = "") else "nlminb\n")
    cat("Convergence:", x$conv, "\n")
    cat("max(|grad|):", format.pval(x$max.sc, digits = 2, eps = 1e-06), "\n")
    cat("\n\n")
    invisible(x)
}
