summary.mnl <- function (object,...){
	b <- coef(object)
	std.err <- sqrt(diag(solve(-object$hessian)))
  	z <- b/std.err
  	p <- 2*(1-pnorm(abs(z)))
  	CoefTable <- cbind(b,std.err,z,p)
  	colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  	lratio <- 1-(object$logLik/object$logLik0)
  	n <- nrow(object$fitted.values)
	out <- list(logLik=object$logLik, logLik0=object$logLik0,
				lratio=lratio, aic=object$aic, bic=object$bic,
				n=n, iter=object$iter[1],
				c=object$convergence,
				message=object$message,
				method=object$method, time=object$time,
				suggested.t = sqrt(log(n)),
				coefs=CoefTable,
				mnl.spec=object$specification
			)
	class(out) <- 'summary.mnl'
	return(out)
}

print.summary.mnl <- function (x, digits = max(3, getOption("digits") - 3), ...){
  cat('\n')
  cat(x$message,'\n\n')
  cat("Log-Likelihood:\t\t\t",signif(x$logLik, digits),"\n")
  cat("Null Log-Likelihood:\t",signif(x$logLik0, digits),"\n")
  cat("Likelihood ratio index:\t", signif(x$lratio, digits), '\n')
  cat('AIC:\t\t\t\t\t', x$aic, '\n')
  cat('BIC:\t\t\t\t\t', x$bic, '\n')
  cat('Sample size:\t\t\t', x$n, '\n')
  cat('Iterations:\t\t\t\t', x$iter, '\n')
  cat('Suggested |t-value| > \t', x$suggested.t, '\n')
  cat('Convergence statistics:\t', x$c, '\n')
  cat("\nEstimated using", x$method, "in", x$time, "s.\n")
  cat("\nCoefficients :\n")
  print(x$coef, digits=digits)
  
  summary(x$mnl.spec)
}
