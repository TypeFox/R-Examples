rmvnormpost <-
function (n, prior.mean, prior.precision, sample.mean, sample.precision, RNORM.METHOD)
{
	post.cov = solve (as.matrix(prior.precision + n*sample.precision))
	post.mean = as.vector (post.cov %*% (prior.precision %*% prior.mean + n*sample.precision %*% sample.mean))
	
	result = rMVNorm (1, mean=post.mean, sigma=post.cov, method=RNORM.METHOD)
#	result = rmvnorm (1, mean=post.mean, sigma=post.cov, method=RNORM.METHOD)
	
	return (result)
}

