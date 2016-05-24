# The 'starma' function estimates the parameters of a STARMA model given
# a set of space time observations, using Kalman filters encoded in C++.
# It also computes the residuals, the log likelihood and the BIC of the model.

# Args:
# -	data: A matrix or a data frame with the observations of the space time
#		process.
# -	wlist: The list of lagged weight matrices, the first one being the
#		 identity.
# -	ar: Either an integer specifying the AR maximum time lag, or
#	    a 1/0 matrix determining if 'row'-th tlag, 'col'-th slag AR 
#	    parameter should be estimated.
# -	ma: Either an integer specifying the MA maximum time lag, or
#	    a 1/0 matrix determining if 'row'-th tlag, 'col'-th slag MA
#	    parameter should be estimated

# Returns:
# A list containing:
#	$phi: The estimated AR parameters
#	$phi_sd: The corresponding standard errors
#	$theta: The estimated MA parameters
#	$theta_sd: The corresponding standard errors
#	$sigma2: The estimated white noise variance matrix
#		   Note that, to achieve parcimony, only the mean of the diagonal
#		   elements should be kept (since the noise is supposed to be 
#		   Gaussian anyway)
#	$residuals: The estimated residuals of the model
#	$loglik: The conditional log likelihood of the model
#	$bic: The corresponding BIC
#	$call: The function call
#	$df: Degrees of freedom of the model: (nb of obs) - (nb of parameters)

starma <- function(data, wlist, ar, ma, iterate=1) UseMethod("starma")

starma.default <- function(data, wlist, ar, ma, iterate=1) {

	if (is.data.frame(data))
		data <- as.matrix(data)

	if (!is.matrix(ar))
		ar <- matrix(1, ar, length(wlist))
	
	if (!is.matrix(ma))
		ma <- matrix(1, ma, length(wlist))

	model <- starmaCPP(data, wlist, ar, ma, iterate)

	# Name the parameter matrices
	if (nrow(model$phi)) {
		model$phi <- replace(model$phi, which(model$phi == 0), NA)
		colnames(model$phi) <- paste("slag", 1:ncol(model$phi) - 1)
		rownames(model$phi) <- paste("tlag", 1:nrow(model$phi))
		model$phi_sd <- replace(model$phi_sd, which(model$phi_sd == 0), NA)
		colnames(model$phi_sd) <- paste("slag", 1:ncol(model$phi_sd) - 1)
		rownames(model$phi_sd) <- paste("tlag", 1:nrow(model$phi_sd))
	}

	if (nrow(model$theta)) {
		model$theta <- replace(model$theta, which(model$theta == 0), NA)
		colnames(model$theta) <- paste("slag", 1:ncol(model$theta) - 1)
		rownames(model$theta) <- paste("tlag", 1:nrow(model$theta))
		model$theta_sd <- replace(model$theta_sd, which(model$theta_sd == 0), NA)
		colnames(model$theta_sd) <- paste("slag", 1:ncol(model$theta_sd) - 1)
		rownames(model$theta_sd) <- paste("tlag", 1:nrow(model$theta_sd))
	}
	
	# As a 'starma' class object
	model$call <- match.call()
	model$df <- nrow(data) * ncol(data) - sum(ar) - sum(ma)
	class(model) <- "starma"

	# Returns the model
	model

}

print.starma <- function(x, ...) {

	cat("Call:\n")
	print(x$call)

	cat("---------------------------------------------------------------\n")

	if (nrow(x$phi)) {
		cat("\nAR parameters:\n")
		printCoefmat(x$phi)
		cat("\nStandard deviation:\n")
		printCoefmat(x$phi_sd)

	cat("---------------------------------------------------------------\n")
	}
	
	if (nrow(x$theta)) {
		cat("\nMA parameters:\n")
		printCoefmat(x$theta)
		cat("\nStandard deviation:\n")
		printCoefmat(x$theta_sd)

	cat("---------------------------------------------------------------\n")
	}

	cat("\nsigma2 estimated as ", 
		sum(diag(x$sigma2))/ncol(x$sigma2), 
		":  log likelihood = ", x$loglik, "  bic = ", x$bic, 
		"\n", sep="")

}
