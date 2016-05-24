#Simulate samples from a multinomial distribution
#The algorithm is to treat each variable as if it
#were binomial at a time

rmultin<-function(n, p)
{
	l<-length(p)
	x<-rep(0,l)
	
	pdenom<-1
	m<-n

	for(i in 1:(l-1) )
	{
		pr<-p[i]/pdenom
		x[i]<-rbinom(1, m, pr)

		m<-m-x[i]
		pdenom<-pdenom-p[i]
	}
	x[l]<-m

	return(x)
}
