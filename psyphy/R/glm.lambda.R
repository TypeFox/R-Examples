glm.lambda <- function(formula, data, NumAlt = 2,
			lambda = seq(0, 0.1, len = 40), 
			plot.it = FALSE, ...) {
	if (missing(data)) 
        data <- environment(formula)
	lam.lst <- vector("list", length(lambda))
	for (ix in seq(length(lambda))) {
		lam.lst[[ix]] <- glm(formula, family = 
			binomial(probit.lambda(NumAlt, lambda[ix])),
			data = data, ...)
		}
	dev <- unlist(sapply(lam.lst, deviance))
	n <- 5
	mdev <- which.min(dev)
	if (mdev > n) {
		ss <- seq(-n, n) + mdev } else {
		ss <- seq(1, 2*n+1) }
	lam.quad <- lm(dev[ss] ~ lambda[ss] + I(lambda[ss]^2))
	lmin <- -0.5 * coef(lam.quad)[2]/coef(lam.quad)[3]
	res <- glm(formula, family = 
			binomial(probit.lambda(NumAlt, lmin)),
			data = data, ...)
	class(res) <- c("lambda", "glm", "lm")
	res$lambda <- as.vector(lmin)
	res$df.residual <- res$df.residual - 1
	res$profile <- data.frame(lambda = lambda, deviance = dev)
	if (plot.it) {
		plot(lambda, dev, xlab = expression(lambda), cex = 0.2,
				ylab = "Deviance", cex.lab = 1.35, type = "p")
		if (length(lambda[ss]) == length(predict(lam.quad))) {
			lines(lambda[ss], predict(lam.quad), lwd = 2) 
#		lines(lambda, predict(lam.quad), lwd = 2)
		abline(v = lmin, lty = 3) }
		}
	cat("\n", "lambda = ", lmin, "\n")	
	res
}	

summary.lambda <- function(object, ...) {
	ans <- summary.glm(object, ...)	
	ans$lambda <- object$lambda
	if (!is.null(object$gam)) ans$gamma <- object$gam
	if(!is.null(object$SEgam)) ans$SEgam <- object$SEgam
	if(!is.null(object$SElambda)) ans$SElambda <- object$SElambda
	class(ans) <- c("summary.lambda")
	return(ans)
}

print.summary.lambda <- function(x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor, 
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
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", 
        format(x$dispersion), ")\n\n", apply(cbind(paste(format(c("Null", 
            "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance", 
            "deviance")]), digits = max(5, digits + 1)), " on", 
            format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"), 
            1, paste, collapse = " "), sep = "")
    if (nchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4, digits + 1)), 
        "\n\n", "Number of Fisher Scoring iterations: ", x$iter, 
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
    
    cat("lambda\t", round(x$lambda, digits), "\t")
    if (!is.null(x$gam)) cat("gamma\t", round(x$gam, digits), "\n") else
    	cat("\n")
	if (!is.null(x$SEgam)) {
		cat("+/-SE(lambda) = \t(", plogis(qlogis(x$lambda) + c(-x$SElambda, x$SElambda)), ")\n")
		cat("+/-SE(gamma) = \t(", plogis(qlogis(x$gam) + c(-x$SEgam, x$SEgam)), ")\n")	
	}
    invisible(x)
	
}
