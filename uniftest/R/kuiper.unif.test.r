kuiper.unif.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
v<-max(x-c(0:(n-1))/n)+max(c(1:n)/n-x)
for(i in 1:nrepl)
{
	z<-runif(n)
	V<-max(z-c(0:(n-1))/n)+max(c(1:n)/n-z)
	if (V>v) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(V=v), p.value=p.value, method="Kuiper test for uniformity",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}