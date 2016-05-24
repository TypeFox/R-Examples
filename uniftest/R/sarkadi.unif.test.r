sarkadi.unif.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
d<-(x-c(1:n)/(n+1))/(c(1:n)*(n-c(1:n)+1))
j<-n*n*sum(d^2)-n*(sum(d))^2
for(i in 1:nrepl)
	{
	z<-runif(n)
	z<-sort(z)
	D<-(z-c(1:n)/(n+1))/(c(1:n)*(n-c(1:n)+1))
	J<-n*n*sum(d^2)-n*(sum(d))^2
	if (J>j) l=l+1
	}
p.value<-l/nrepl
RVAL<-list(statistic=c(J=j), p.value=p.value, method="Sarkadi-Kosik test for uniformity",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}