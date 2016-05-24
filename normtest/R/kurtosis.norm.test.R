kurtosis.norm.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x1<-sum(x)/n
t<-n*(sum((x-x1)^4))/((sum((x-x1)^2))^2)
for(i in 1:nrepl)
{
	z<-rnorm(n)
	z1<-sum(z)/n
	T<-n*(sum((z-z1)^4))/((sum((z-z1)^2))^2)
	if (abs(T-3)>abs(t-3)) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(T=t), p.value=p.value, method="Kurtosis test for normality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}