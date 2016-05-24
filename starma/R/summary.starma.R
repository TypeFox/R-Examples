# The 'summary' method for the 'starma' class outputs the estimated 
# coefficients of a model and the corresponding p-values.

summary.starma <- function(object, ...) {

	# Get coefficient values and compute p-values
	coefs <- c(as.vector(t(object$phi)), as.vector(t(object$theta)))
	sd <- c(as.vector(t(object$phi_sd)), as.vector(t(object$theta_sd)))
	tval <- coefs / sd
	pval <- 2 * pt(-abs(tval), df = object$df)

	# Form the aggregated matrix of coefficients
	tab <- data.frame(Estimate = coefs,
			 	"Std. Error" = sd,
			 	t.value = tval,
			 	p.value = pval)

	# Pretty-labels the rows
	label.phi <- NULL
	label.theta <- NULL
	if (nrow(object$phi)) {
		tlag.phi <- matrix(1:nrow(object$phi), 
					 ncol(object$phi), nrow(object$phi), byrow=T)
		slag.phi <- matrix(1:ncol(object$phi) - 1, 
					 ncol(object$phi), nrow(object$phi))
		label.phi <- paste("phi", tlag.phi, slag.phi, sep="")
	}
	if (nrow(object$theta)) {
		tlag.theta <- matrix(1:nrow(object$theta), 
					  ncol(object$theta), nrow(object$theta), byrow=T)
		slag.theta <- matrix(1:ncol(object$theta) - 1, 
					   ncol(object$theta), nrow(object$theta))
		label.theta <- paste("theta", tlag.theta, slag.theta, sep="")
	}
	rownames(tab) <- c(label.phi, label.theta)

	# Remove NA rows
	tab <- tab[complete.cases(tab), ]
	
	# Outputs 'summary.starma' class object
	out <- list(call=object$call,
			coefficients=tab)

	class(out) <- "summary.starma"
	out

}

print.summary.starma <- function(x, ...) {

	cat("Call:\n")
	print(x$call)
	cat("\n")
	
	printCoefmat(x$coefficients, P.values=T, has.Pvalue=T)

}