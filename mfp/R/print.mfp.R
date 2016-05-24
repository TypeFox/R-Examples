print.mfp <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }
    if (!is.null(x$fail)) {
        cat(" MFP failed.", x$fail, "\n")
        return()
    }
#
	cox <- (x$family$family=="Cox")
#
# Deviance table
#
		cat("\nDeviance table:")
		cat("\n \t\t Resid. Dev")
		cat("\nNull model\t", x$dev.null)
		cat("\nLinear model\t", x$dev.lin)
		cat("\nFinal model\t", x$dev)
		cat("\n")
#
# Fractional Polynomials
#
    cat("\nFractional polynomials:\n")
    fptable(x)
#
    cat("\n")
#
# Transformations used
#
	cat("\nTransformations of covariates:\n"); print(x$trafo); cat("\n")

#
    savedig <- options(digits = digits)
    on.exit(options(savedig))
if(cox) {
    coef <- x$coef
    se <- sqrt(diag(x$var))
    if(!cox) if (is.null(coef) | is.null(se)) 
        stop("Input is not valid")
    if(!is.null(coef)) {
	 if(x$rescale & any(x$scale[,1]>0)) {
			cat("Re-Scaling:\nNon-positive values in some of the covariates. No re-scaling was performed.\n\n")
			x$rescale <- FALSE
		}
#
	 if (is.null(x$naive.var)) {
        tmp <- cbind(coef, exp(coef), se, coef/se, signif(1 - 
            pchisq((coef/se)^2, 1), digits - 1))
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
            "se(coef)", "z", "p"))
     }
     else {
        nse <- sqrt(diag(x$naive.var))
        tmp <- cbind(coef, exp(coef), nse, se, coef/se, signif(1 - 
            pchisq((coef/se)^2, 1), digits - 1))
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
            "se(coef)", "robust se", "z", "p"))
     }
    prmatrix(tmp)
    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    if (is.null(x$df)) 
        df <- sum(!is.na(coef))
    else df <- round(sum(x$df), 2)
    cat("\n")
    cat("Likelihood ratio test=", format(round(logtest, 2)), 
        "  on ", df, " df,", " p=", format(1 - pchisq(logtest, 
            df)), sep = "")
	}
	if(is.null(coef)) {
    cat("\n")
    cat("Null model\n log likelihood=", format(round(x$loglik, 2)), "\n")
    }	
    omit <- x$na.action
    if (length(omit)) 
        cat(" n=", x$n, " (", naprint(omit), ")\n", sep = "")
    else cat(" n=", x$n, "\n")
    if (length(x$icc)) 
        cat("   number of clusters=", x$icc[1], "    ICC=", format(x$icc[2:3]), 
            "\n")
    invisible(x)
} else {
    if (length(coef(x))) {
		if(x$rescale & any(x$scale[,1]>0)) {
			cat("Re-Scaling:\nNon-positive values in some of the covariates. No re-scaling was performed.\n\n")
			x$rescale <- FALSE
		}
        if(x$rescale) cat("Rescaled coefficients") else cat("Coefficients")
        if (is.character(co <- x$contrasts)) 
            cat("  [contrasts: ", apply(cbind(names(co), co), 
                1, paste, collapse = "="), "]")
        cat(":\n")
        print.default(format(x$coefficients, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", 
        x$df.residual, "Residual\n")
    cat("Null Deviance:\t   ", format(signif(x$null.deviance, 
        digits)), "\nResidual Deviance:", format(signif(x$deviance, 
        digits)), "\tAIC:", format(signif(x$aic, digits)), "\n")
    invisible(x)
}
}


fptable <- function(x) print(x$fptable)
