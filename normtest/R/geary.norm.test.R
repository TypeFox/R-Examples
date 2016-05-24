geary.norm.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x1<-sum(x)/n
s<-sqrt(sum((x-x1)^2))
d<-sum(abs(x-x1))/(s*sqrt(n))
for(i in 1:nrepl)
{
	z<-rnorm(n)
	z1<-sum(z)/n
	s<-sqrt(sum((z-z1)^2))
	D<-sum(abs(z-z1))/(s*sqrt(n))
	if (D>d) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(d=d), p.value=p.value, method="Geary test for normality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}