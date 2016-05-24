frosini.norm.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
y<-mean(x)
s<-sqrt(var(x)*(n-1)/n)
b<-(1/sqrt(n))*sum(abs(pnorm((x-y)/s)-(c(1:n)-0.5)/n))
for(i in 1:nrepl)
{
	z<-rnorm(n)
	z<-sort(z)
	y<-mean(z)
	s<-sqrt(var(z)*(n-1)/n)
	B<-(1/sqrt(n))*sum(abs(pnorm((z-y)/s)-(c(1:n)-0.5)/n))
	if (B>b) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(B=b), p.value=p.value, method="Frosini test for normality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}