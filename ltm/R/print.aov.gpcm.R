print.aov.gpcm <-
function (x, digits = 3, ...) {
    if (!inherits(x, "aov.gpcm"))
        stop("Use only with 'aov.gpcm' objects.\n")
    p.val <- round(x$p.value, 3)
    p.val <- if (p.val < 0.001) "<0.001" else p.val
    dat <- data.frame(AIC = round(c(x$aic0, x$aic1), 2), BIC = round(c(x$bic0, x$bic1), 2), 
                log.Lik = round(c(x$L0, x$L1), 2), LRT = c(" ", round(x$LRT, 2)), df = c(x$nb0, x$nb1), 
                p.value = c("", p.val), row.names = c(x$nam0, x$nam1))
    if (x$simulate.p.value)
        cat("\n Likelihood Ratio Table (based on parametric Bootstrap with", x$B,  "datasets)\n\n")
    else
        cat("\n Likelihood Ratio Table\n\n")
    print(dat)
    cat("\n\n")
    invisible(x)
}
