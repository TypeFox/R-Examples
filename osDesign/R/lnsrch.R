lnsrch <-
function(gamm0,f0,g0,p,nstrata,nn1,nn0,n1a,n0a,grpa,repp,xx,ofs,yy)
{
	scale <- 1.0
	rel.error <- max(abs(p)/pmax(gamm0,0.1))
	if(rel.error < 1e-07)
		gamm <- gamm0
	else
	{
		gamm <- gamm0 + p
		l <- lagrange(gamm,nstrata,nn1,nn0,n1a,n0a,grpa,repp,xx,ofs,yy)
		f <- 0.5*sum(l^2)
		slope <- sum(g0*p)
		# backtrack if the function does not decrease sufficiently
		if(f > f0+0.0001*slope)
		{
			scale <- -slope/(2*(f-f0-slope))
			scale <- min(scale,0.5)
			scale <- max(0.05,scale)
    }
		gamm <- gamm0 + scale*p
	}
	gamm
}
