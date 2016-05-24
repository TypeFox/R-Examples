print.thregIcure <-
 function(x, digits=max(options()$digits - 4, 3), ...)
 {
        if (!is.null(cl<- x$call)) {
		cat("Call:\n")
		dput(cl)
		cat("\n")
	}

        coef <- x$coefficients
        se <- sqrt(diag(x$var))
	tmp <- cbind(coef, se, coef/se,
	       signif(1 - pchisq((coef/ se)^2, 1), digits -1))
	dimnames(tmp) <- list(names(coef), c("coef",
	    "se(coef)", "z", "p"))
        prmatrix(tmp)
	cat("\n")
	cat("Log likelihood =", format(round(x$loglik, 2)), ",", " AIC =", format(round(x$AIC, 2)), ", ", "Goodness of fit test: p-value =", format(round(x$Pvalue, 2)), sep="")
	cat("\n")
        invisible(x)

}


