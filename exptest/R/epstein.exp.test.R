epstein.exp.test<-function(x, simulate.p.value=FALSE, nrepl=2000) 
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
x<-c(0,x)
D<-(n:1)*(x[2:(n+1)]-x[1:n])
t<-2*n*(log(sum(D)/n)-(sum(log(D)))/n)/(1+(n+1)/(6*n))
l<-0
if (simulate.p.value)
{
	for(i in 1:nrepl)
	{
	z<-rexp(n)
	z<-sort(z)
	z<-c(0,z)
	D<-(n:1)*(z[2:(n+1)]-z[1:n])
	T<-2*n*(log(sum(D)/n)-(sum(log(D)))/n)/(1+(n+1)/(6*n))
	if (T>t) l=l+1
	}
	p.value<-l/nrepl
}
else
{
	p.value<-1-pchisq(t,n-1)
}
RVAL<-list(statistic=c(EPS=t), p.value=p.value, method="Epstein test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}