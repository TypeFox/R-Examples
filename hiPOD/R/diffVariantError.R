diffVariantError <-
function(Xmean, N.p, error, N.test=1)
{
	theta <- (1-error)/(2*N.p) + error*(1-1/(2*N.p))

	for(v in 1:Xmean)
	{
		p1 <- pbinom(v-1, Xmean, theta, lower.tail=FALSE)
		p0 <- pbinom(v-1, Xmean, error, lower.tail=FALSE)
		
		if(p1/p0 > 5 & p0<0.05/N.test) break;
	}
	
	c(v, p0, p1)
}

