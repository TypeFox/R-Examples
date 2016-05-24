print.aov.tpm <-
function (x, digits = 3, ...) {
    if (!inherits(x, "aov.tpm"))
        stop("Use only with 'aov.tpm' objects.\n")
    p.val <- round(x$p.value, 3)
    p.val <- if (p.val < 0.001) "<0.001" else p.val
    dat <- data.frame(AIC = round(c(x$aic0, x$aic1), 2), BIC = round(c(x$bic0, x$bic1), 2), 
                log.Lik = round(c(x$L0, x$L1), 2), LRT = c(" ", round(x$LRT, 2)), df = c("", x$df), 
                p.value = c("", p.val), row.names = c(x$nam0, x$nam1))
    cat("\n Likelihood Ratio Table\n")
    print(dat)
    cat("\n\n")
    invisible(x)
}
