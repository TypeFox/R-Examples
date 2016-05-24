spiegelhalter.norm.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
x1<-sum(x)/n
s<-sqrt(((sum((x-x1)^2))/(n-1)))
u<-(x[n]-x[1])/s
c<-factorial(n)^(1/(n-1))/(2*n)
g<-(sum(abs(x-x1)))/(s*sqrt(n*(n-1)))
t<-((c*u)^(1-n)+g^(1-n))^(1/(n-1))
for(i in 1:nrepl)
{
	z<-rnorm(n)
	z<-sort(z)
	z1<-sum(z)/n
	s<-sqrt((sum((z-z1)^2))/(n-1))
	u<-(z[n]-z[1])/s
	g<-sum(abs(z-z1))/(s*sqrt(n*(n-1)))
	T<-((c*u)^(1-n)+g^(1-n))^(1/(n-1))
	if (T>t) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(T=t), p.value=p.value, method="Spiegelhalter test for normality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}