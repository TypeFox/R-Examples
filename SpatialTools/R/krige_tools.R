krige.uk <- function(y, V, Vp, Vop, X, Xp, nsim = 0, Ve.diag = NULL, method = "eigen")
{
  # Check arguments.
	ins <- krige_arg_check(y, V, Vp, Vop, X = NULL, Xp = NULL, m = 0, nsim, Ve.diag, method)

	###compute matrix products for future use
	ViX <- solve(V, X)
	XtViX <- crossprod(ViX, X)

	#compute gls estimates of regression coefficients
	coeff <- solve(XtViX, crossprod(ViX, y))

	#compute kriging weights
	w <- solve(V, Vop - X %*% solve(XtViX, crossprod(ViX, Vop) - t(Xp)))
	
	#blup for yp
	pred <- crossprod(w, y)

	#variance of (yp - pred)
	mspe <- colSums((V %*% w) * w) - 2 * colSums(w * Vop) + diag(Vp)

	out <- list(pred = pred, mspe = mspe, coeff = coeff,
	vcov.coef = solve(XtViX))

	# generate conditional realizations if nsim > 0
	if(nsim > 0)
	{
		# determine number of observed and predicted data values
		n <- nrow(V); np <- nrow(Vp)

		# modify V to account for error
		Vomod <- V - diag(Ve.diag)

		# Combine observed and predicted covariance matrices
		Va <- rbind(cbind(Vomod, Vop), cbind(t(Vop), Vp))

		# Generate mean zero normal random variables with appropriate covariance
		# matrix (does not include error)
		newsim <- decomp.cov(Va, method = method) %*%
		matrix(stats::rnorm(nrow(Va) * nsim), ncol = nsim)

		# Take difference of observed data and simulated observations
		# (now including error)
		newZ <- matrix(y, nrow = n, ncol = nsim) - newsim[1:n,] +
		matrix(stats::rnorm(n * nsim, sd = sqrt(Ve.diag)), nrow = n, ncol = nsim)

		# Create conditional realizations
		sim <- newsim[-(1:n),] + crossprod(w, newZ)

		out$sim <- sim
		class(out) <- "krigeConditionalSample"
	}
	return(out)
}


krige.ok <- function(y, V, Vp, Vop, nsim = 0, Ve.diag = NULL, method = "eigen")
{
	# Check arguments.
	ins <- krige_arg_check(y, V, Vp, Vop, X = NULL, Xp = NULL, m = 0, nsim, Ve.diag, method)

	out <- .Call( "krige_ok", ys = y, Vs = V, Vps = Vp, Vops = Vop, 
		nsims = nsim, Vediags = ins$Ve.diag, method = ins$method, 
		PACKAGE = "SpatialTools")

	#convert one-dimensional matrices to vectors
	out$pred <- as.vector(out$pred)	
	out$mspe <- as.vector(out$mspe)
	out$coeff <- as.vector(out$coeff)
	if(nsim > 0) 		class(out) <- "krigeConditionalSample"
	return(out)
}


krige.sk <- function(y, V, Vp, Vop, m = 0, nsim = 0, Ve.diag = NULL, method = "eigen")
{
	ins <- krige_arg_check(y, V, Vp, Vop, X = NULL, Xp = NULL, m = m, nsim, Ve.diag, method)

	out <- .Call( "krige_sk", ys = y, Vs = V, Vps = Vp, Vops = Vop, ms = m, 
		nsims = nsim, Vediags = ins$Ve.diag, method = ins$method, 
		PACKAGE = "SpatialTools")
	
	#convert one-dimensional matrices to vectors
	out$pred <- as.vector(out$pred)	
	out$mspe <- as.vector(out$mspe)
	if(nsim > 0) 		class(out) <- "krigeConditionalSample"
	return(out)
}