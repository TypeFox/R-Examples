skewness.norm.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x1<-sum(x)/n
t<-sqrt(n)*(sum((x-x1)^3))/((sum((x-x1)^2))^(3/2))
for(i in 1:nrepl)
{
	z<-rnorm(n)
	z1<-sum(z)/n
	T<-sqrt(n)*(sum((z-z1)^3))/((sum((z-z1)^2))^(3/2))
	if (abs(T)>abs(t)) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(T=t), p.value=p.value, method="Skewness test for normality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}