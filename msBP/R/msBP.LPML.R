msBP.LPML <- function(G0_xi, g0_xi, postW, base=FALSE, mu=0, sigma2=1, x=NULL)
{
	
	n <- length(G0_xi)
	R <- nrow(postW)
	lik <- matrix(NA, R, n)
	if(base)
	{
		for(i in 1:R)
		{
			W <- vec2tree(postW[i,])
			lik[i, ] <- dnorm(x, mu[i], sqrt(sigma2[i])) * msBP.pdf(W, n, pnorm(x, mu[i], sqrt(sigma2[i])))$dens
		}
	
	}
	else
	{
		for(i in 1:R)
		{
			W <- vec2tree(postW[i,])
			lik[i, ] <- g0_xi * msBP.pdf(W, n, G0_xi)$dens
		}
	}	
	LPML <- sum(log(1/apply(1/lik, 2, mean)))		
	meanCPO <- mean(log(1/apply(1/lik, 2, mean)))
	medianCPO <- median(log(1/apply(1/lik, 2, mean)))
	list(LPML=LPML, elogCPO=meanCPO, mlogCPO=medianCPO)
}