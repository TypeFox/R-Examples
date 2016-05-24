print.mrqap.dsp <- function(x, ...) {
    cat("MRQAP with Double-Semi-Partialing (DSP)\n\n")
    cat(paste("Formula: ", Reduce(paste, deparse(x$formula)), "\n\n"))
    cat("Coefficients:\n")
    results_mat <- as.vector(format(as.numeric(x$coefficients)))
    results_mat <- cbind(results_mat, as.vector(format(x$P.lesser)))
    results_mat <- cbind(results_mat, as.vector(format(x$P.greater)))
    results_mat <- cbind(results_mat, as.vector(format(x$P.values)))
    colnames(results_mat) <- c("Estimate", "P(\u03B2>=r)", "P(\u03B2<=r)", "P(|\u03B2|<=|r|)")
    rownames(results_mat) <- as.vector(names(x$P.values))
    print.table(results_mat)
    if (names(x$P.values)[1] == "intercept") { 
    	mss <- sum((fitted(x) - mean(fitted(x)))^2)
    	df.int <- 1
    } else {
    	mss <- sum(fitted(x)^2)
    	df.int <- 0
    }
    rss <- sum(resid(x)^2)
    resvar <- rss/x$df.residual
    f_stat <- c(value = (mss/(x$rank - df.int))/resvar, numdf = x$rank - df.int, dendf = x$df.residual)
    r.squared <- mss/(mss + rss)
    adj.r.squared <- 1 - (1 - r.squared) * ((x$n - df.int)/x$df.residual)
    sigma <- sqrt(resvar)
    cat("\nResidual standard error:", format(sigma, digits = 4), "on", x$df.residual, "degrees of freedom\n")
    cat("F-statistic:", formatC(f_stat[1], digits = 4), "on", f_stat[2], "and", f_stat[3], "degrees of freedom, p-value:", formatC(1 - pf(f_stat[1], f_stat[2], f_stat[3]), digits = 4), "\n")
    cat("Multiple R-squared:", format(r.squared, digits = 4), "\t")
    cat("Adjusted R-squared:", format(adj.r.squared, digits = 4), "\n")
    cat("AIC:", x$AIC, "\n")
    cat("\n")
}
