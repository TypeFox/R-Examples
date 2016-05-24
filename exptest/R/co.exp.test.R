co.exp.test<-function(x, simulate.p.value=FALSE, nrepl=2000)
{
DNAME <- deparse(substitute(x))
n<-length(x)
y<-x/mean(x)
y<-log(y)*(1-y)
co<-sum(y)+n
v<-sqrt(6/n)*co/pi
l<-0
if (simulate.p.value)
{
	for(i in 1:nrepl)
	{
	z<-rexp(n)
	z<-z/mean(z)
	z<-log(z)*(1-z)
	CO<-sum(z)+n
	V<-sqrt(6/n)*CO/pi
	if (abs(V)>abs(v)) l=l+1
	}
	p.value<-l/nrepl
}
else
{
	p.value<-2*(1-pnorm(abs(v)))
}
RVAL<-list(statistic=c(COn=co), p.value=p.value, method="Test for exponentiality based on the statistic of Cox and Oakes",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}