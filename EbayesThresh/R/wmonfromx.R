"wmonfromx" <-
function(xd, prior = "laplace", a = 0.5, tol = 1e-008, maxits = 20)
{
#
#   Find the monotone marginal maximum likelihood estimate of the mixing weights
#    for the Laplace prior with parameter a.  It is assumed that the 
#    noise variance is equal to one.
#
#  Find the beta values and the minimum weight
#
	pr <- substring(prior, 1, 1)
	nx <- length(xd)
	wmin <- wfromt(sqrt(2 * log(length(xd))), prior, a)
	winit <- 1
	if(pr == "l")
		beta <- beta.laplace(xd, a)
	if(pr == "c") beta <- beta.cauchy(xd)	#
#   now conduct iterated weighted least squares isotone regression
#    
	w <- rep(winit, length(beta))
	for(j in (1:maxits)) {
		aa <- w + 1/beta
		ps <- w + aa
		ww <- 1/aa^2
		wnew <- isotone(ps, ww, increasing = FALSE)
		wnew <- pmax(wmin, wnew)
		wnew <- pmin(1, wnew)
		zinc <- max(abs(range(wnew - w)))
		w <- wnew
		if(zinc < tol)
			return(w)
	}
	cat("Warning: more iterations required to achieve convergence \n")
	return(w)
}
