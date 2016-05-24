permut<-function(n)
{
	ind<-1:n
	a<-sample(ind,n)
	asort <- sort(a)
	b <- sample(asort,n)
	ind <- b
 	return(ind)
}
