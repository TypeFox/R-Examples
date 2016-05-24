llr<-function(x,p,q,r,m, lawley.correction=F)
{
	N <- sum(x)	
	ll <- 2*( x*log(r/p)+x*log(r/q) )
	if (lawley.correction)
	{
		ll <- ll / (1+ (sum(N/x)-1) /(3*N*(factorial(m)-1)))
	}
	
	return( ll );
	
}