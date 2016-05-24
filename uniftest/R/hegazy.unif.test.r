hegazy.unif.test<-function(x, nrepl=2000, p=1)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
t<-sum((abs(x-c(1:n)/(n+1)))^p)/n
for(i in 1:nrepl)
{
	z<-runif(n)
	z<-sort(z)
	T<-sum((abs(z-c(1:n)/(n+1)))^p)/n
	if (T>t) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(T=t), p.value=p.value, method="Hegazy-Green test for uniformity",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}