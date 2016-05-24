`print.summary.cusp` <-
function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("Deviance Residuals: \n")
    if (x$df.residual > 5) {
        x$deviance.resid <- quantile(x$deviance.resid, na.rm = TRUE)
        names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", 
            "Max")
    }
    print.default(x$deviance.resid, digits = digits, na.print = "", 
        print.gap = 2)
    if (length(x$aliased) == 0) {
        cat("\nNo Coefficients\n")
    }
    else {
        if (!is.null(df <- x$df) && (nsingular <- df[3] - df[1])) 
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
                sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if (!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    }
    cat("\n\n", apply(cbind(paste(format(c("Null", "Linear", 
        "Logist", x$resid.name), justify = "right"), "deviance:"), 
        format(unlist(x[c("null.deviance", "r2lin.dev", "r2log.dev", 
            "deviance")]), digits = max(5, digits + 1)), " on", 
        format(unlist(x[c("df.null", "r2lin.df", "r2log.df", 
            "df.residual")])), " degrees of freedom\n"), 1, paste, 
        collapse = " "), sep = "")
    inf.crit = data.frame("R Squared" = c(x$r2lin.r.squared, 
        x$r2log.r.squared, x$r2cusp.r.squared), logLik = c(x$r2lin.logLik, 
        x$r2log.logLik, x$r2cusp.logLik), npar = c(x$r2lin.npar, x$r2log.npar,
        x$r2cusp.npar), AIC = c(x$r2lin.aic, 
        x$r2log.aic, x$r2cusp.aic), AICc = c(x$r2lin.aicc, x$r2log.aicc, 
        x$r2cusp.aicc), BIC = c(x$r2lin.bic, x$r2log.bic, x$r2cusp.bic))
    rownames(inf.crit) <- paste(c("Linear", "Logist", "Cusp"), 
        "model")
    cat("\n")
    print(na.omit(inf.crit))
    cat("---\nNote: R.Squared for cusp model is Cobb's pseudo-R^2. This value\n      can become negative.")
    .chi2 = abs(-2 * ( x$r2lin.logLik - x$r2cusp.logLik ))
    .chi2df = abs(x$r2lin.npar - x$r2cusp.npar)
    cat("\n\n\tChi-square test of linear vs. cusp model\n\n")
    cat("X-squared = ", formatC(.chi2), ", df = ", .chi2df, ", p-value = ", formatC(1-pchisq(.chi2, .chi2df)), sep="")
    if (nchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
    cat("\n\n", "Number of optimization iterations: ", x$iter, 
        "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2, 
                  digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\n")
    invisible(x)
}

