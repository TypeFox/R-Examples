print.summary.weibull.frailty <-
function (x, digits = max(4, getOption("digits") - 4), ...) {
    cat("\n\n\tWeibull Relative Risk Model with Gamma Frailty\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat("Data Descriptives:\n")
    cat("Number of groups:", length(unique(x$id)), "\n")
    cat("Number of observations:", length(x$id), "\n")
    cat("Total Number of Events: ", sum(x$d), "\n")
    cat("\nModel Summary:\n")
    model.sum <- data.frame(log.Lik = x$logLik, AIC = x$AIC, BIC = x$BIC, row.names = "")
    print(model.sum)
    if (!is.null(x$coefTab)) {
        cat("\nCoefficients:\n")
        out <- as.data.frame(round(x$coefTab, digits))
        ind <- out$"p-value" == 0
        out$"p-value" <- sprintf(paste("%.", digits, "f", sep = ""), out$"p-value")
        out$"p-value"[ind] <- paste("<0.", paste(rep("0", digits - 1), collapse = ""), "1", sep = "")
        print(out)
    }
    cat("\nShape:", round(x$shape, digits), "\nScale:", round(x$scale, digits), 
        "\nFrailty variance:", round(x$var.frailty, digits))
    cat("\n\n")
    invisible(x)
}
