gini.exp.test<-function(x, simulate.p.value=FALSE, nrepl=2000)
{
DNAME <- deparse(substitute(x))
n<-length(x)
a<-1:(n-1)
b<-(n-1):1
a<-a*b
x<-sort(x)
k<-x[2:n]-x[1:(n-1)]
g<-sum(k*a)/((n-1)*sum(x))
v<-abs(g-0.5)*(sqrt(12*(n-1)))
l<-0
if(simulate.p.value)
{
	for(i in 1:nrepl)
	{
	z<-rexp(n)
	z<-sort(z)
	k<-z[2:n]-z[1:(n-1)]
	G<-sum(k*a)/((n-1)*sum(z))
	V<-abs(G-0.5)*(sqrt(12*(n-1)))
	if (V>v) l=l+1
	}
	p.value<-l/nrepl
}
else
{
	p.value<-2*(1-pnorm(v))
}
RVAL<-list(statistic=c(Gn=g), p.value=p.value, method="Test for exponentiality based on the Gini statistic",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}