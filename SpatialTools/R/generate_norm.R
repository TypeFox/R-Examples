rmvnorm <- function(nsim = 1, mu, V, method = "eigen")
{
	mu <- as.vector(mu)
	
	# check arguments of function
	rmvnorm_arg_check(nsim, mu, V, method)
	
	# decompose covariance matrix
	decomp.V <- decomp.cov(V, method = method)
	
	# return simulated values
	return(mu + decomp.V %*% matrix(stats::rnorm(nrow(V) * nsim), nrow = nrow(V), ncol = nsim))
}

rcondnorm <- function(nsim = 1, y, mu, mup, V, Vp, Vop, method = "eigen")
{
	mc <- mup + crossprod(Vop, solve(V, y - mu))
	Vc <- Vp - crossprod(Vop, solve(V, Vop))
	return(rmvnorm(nsim, mu = mc, V = Vc, method = method))
}

