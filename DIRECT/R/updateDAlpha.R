updateDAlpha <-
function (alpha.curr, prior.shape, prior.rate, alpha.sd, nClust, nSample, method=c("Gibbs", "MH"))
{
	method = match.arg (method)
	if (method=="Gibbs")
	{
		# Sample eta, which has conditional distribution given alpha and n
		# Beta (alpha+1,n)
		eta = rbeta (1, alpha.curr+1, nSample)
		
		# Sample alpha, which has a mixture of two gamma distributions given eta and nClust
		# compute weights
		q = prior.rate - log (eta)
		
		w1 = prior.shape + nClust -1
		w2 = nSample * q
		w = w1+w2
		
		# Draw a sample from the mixture for alpha
		if (runif(1) < (w1/w)) result = rgamma (1, shape=prior.shape+nClust, rate=q)
		else result = rgamma (1, shape=prior.shape+nClust-1, rate=q)
	}
	else
	{
# An MH sampler
		alpha.tmp = rnorm (1, mean=alpha.curr, sd=alpha.sd)
		alpha.new = reflect (alpha.tmp, lower=0, upper=1e8)
		
		target.new = lgamma (alpha.new) - lgamma (alpha.new + nSample) + nClust*log (alpha.new) + dgamma (alpha.new, shape=prior.shape, rate=prior.rate, log=TRUE)
		target.curr = lgamma (alpha.curr) - lgamma (alpha.curr + nSample) + nClust*log (alpha.curr) + dgamma (alpha.curr, shape=prior.shape, rate=prior.rate, log=TRUE)
		
		accept.prob.log = min (0, target.new - target.curr)
		result = ifelse (log (runif(1)) < accept.prob.log, alpha.new, alpha.curr)
	}
	
	return (result)
}

