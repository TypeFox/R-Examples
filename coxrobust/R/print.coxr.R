
print.coxr <- function(x, digits = max(3, getOption("digits") - 4), ...) {

    if ( !inherits(x, "coxr") ) {
        stop("use only with \"coxr\" objects")
    }

    op <- options(digits=digits)
    on.exit(options(op))

    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")

    sd.ple <- sqrt( diag( x$var.ple ) )
    sd <- sqrt(diag(x$var))
    df <- sum(!is.na(x$coef))

    tmp <- cbind(x$ple.coefficients, exp(x$ple.coefficients), sd.ple,
                2*(1-pnorm(abs(x$ple.coefficients/sd.ple))))
    dimnames(tmp)[[1]] <- names(x$coef)
    dimnames(tmp)[[2]] <- c("coef", "exp(coef)", "se(coef)", "p")

    cat("Partial likelihood estimator\n")
    print(tmp)
    cat("\nWald test=", x$wald.test, " on ", df, " df,", " p=",
        1 - pchisq(x$wald.test, df), "\n", sep="")

    tmp <- cbind(x$coef, exp(x$coef), sd, 2*(1-pnorm(abs(x$coef/sd))))
    dimnames(tmp)[[1]] <- names(x$coef)
    dimnames(tmp)[[2]] <- c("coef", "exp(coef)", "se(coef)", "p")

    cat("\nRobust estimator\n")
    print(tmp)
    cat("\nExtended Wald test=", x$ewald.test, " on ", df, " df,", " p=",
        1 - pchisq(x$ewald.test, df), "\n\n", sep="")

    invisible(x)

}