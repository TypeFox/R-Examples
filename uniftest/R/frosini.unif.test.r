frosini.unif.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
b<-1/sqrt(n)*sum(abs(x-(c(1:n)-0.5)/n))
for(i in 1:nrepl)
{
	z<-runif(n)
	z<-sort(z)
	B<-1/sqrt(n)*sum(abs(z-(c(1:n)-0.5)/n))
	if (B>b) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(B=b), p.value=p.value, method="Frosini test for uniformity",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}