probDetEqual3 <-
function(MAF, Xmean, T, N.p, error)
{
	pDetEqual <- 0;
	for(n in 0:(2*N.p))
	{
		pConf <- dbinom(n,2*N.p,MAF);
		pDet.not <- 0;
		theta <- (1-error)*n/(2*N.p) + error*(1-n/(2*N.p))
		
		pDet <- pbinom(T-1, Xmean, theta, lower.tail=FALSE)

		pDetEqual <- pDetEqual + pConf*pDet;
	}	
	pDetEqual
}

