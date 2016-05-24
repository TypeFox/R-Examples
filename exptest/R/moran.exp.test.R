moran.exp.test<-function(x, simulate.p.value=FALSE, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
y<-mean(x)
t<--digamma(1) + mean(log(x/y))
if (simulate.p.value)
{
	for(i in 1:nrepl)
	{
	z<-rexp(n)
	y<-mean(z)
	t<--digamma(1) + mean(log(z/y))
	if (abs(T)>abs(t)) l=l+1
	}
	p.value<-l/nrepl
}
else
{
	p.value<-2*(1-pnorm(abs(sqrt(n)*t/sqrt(((pi^2)/6)-1))))
}
RVAL<-list(statistic=c(T=t), p.value=p.value, method="Moran test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}