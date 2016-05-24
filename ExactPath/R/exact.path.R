exact.path = function(X, y, max.var=20, verbose=FALSE)
{
	n = nrow(X)
	p = ncol(X)
    X = apply(X, 2, scale)
    X = X/sqrt(n-1)
    y = scale(y)/sqrt(n-1)

	lambda = max(abs(t(X) %*%y)) + 1     # start at a value larger than lambda_max  
	ini.tau = rep(0,p)

    breaks = tau = beta = score = NULL
	count = 0
	while (count < min(max.var, p)+4)
	{
		temp = LASSO.exact(lambda, ini.tau, X, y)    
		ttt = abs(ini.tau-temp$tau)>0.5
		if (all(!ttt)) lambda = 0 else lambda= max(temp$breaks[ttt])
		ini.tau = temp$tau
		count = count+1
		if (verbose) {print(temp); cat("\nlambda_",count," = ",lambda,"\n\n\n", sep="")}

		breaks = c(breaks, lambda)
		tau = cbind(tau, temp$tau)
		score = cbind(score, temp$score)
		beta = cbind(beta, temp$beta)
	}
	
	structure(list(breaks=breaks, tau=tau, beta=beta, score=score), class="path")            
}
