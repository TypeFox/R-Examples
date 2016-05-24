hegazy1.norm.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
s<-sqrt(var(x)*(n-1)/n)
x1<-sum(x)/n
x<-(x-x1)/s
b<-qnorm(c(1:n)/(n+1))
t<-sum(abs(x-b))/n
for(i in 1:nrepl)
{
	z<-rnorm(n)
	z<-sort(z)
	s<-sqrt(var(z)*(n-1)/n)
	z1<-sum(z)/n
	z<-(z-z1)/s
	T<-sum(abs(z-b))/n
	if (T>t) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(T=t), p.value=p.value, method="Hegazy-Green test for normality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}