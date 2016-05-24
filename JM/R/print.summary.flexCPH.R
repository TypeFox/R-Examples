print.summary.flexCPH <-
function (x, digits = max(4, getOption("digits") - 4), ...) {
   cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
   cat("Data Descriptives:\n")
   pcEv <- round(100 * sum(x$d) / length(x$d), 1)
   cat("Number of sample units:", length(x$d), "\n")
   cat("Number of Events: ", sum(x$d), " (", pcEv, "%)", sep = "", "\n")
   cat("\nModel Summary:\n")
   cat("B-spline basis internal knots at: ", 
        paste(round(exp(x$knots[-c(1, length(x$knots))]), 2), collapse = ", "), "\n\n", sep = "")
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
   cat("\n\n")
   invisible(x)
}
