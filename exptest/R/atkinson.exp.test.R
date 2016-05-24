atkinson.exp.test<-function(x, p=0.99, simulate.p.value=FALSE, nrepl=2000)
{
DNAME <- deparse(substitute(x))
k<-0
n<-length(x)
y<-mean(x)
m<-mean(x^p)
r<-(m^(1/p))/y
t<-sqrt(n)*abs(r-(gamma(1+p))^(1/p))
if (simulate.p.value)
{
	for(i in 1:nrepl)
	{
	z<-rexp(n)
	y<-mean(z)
	m<-mean(z^p)
	r<-(m^(1/p))/y
	T<-sqrt(n)*abs(r-(gamma(1+p))^(1/p))
	if (T>t) k=k+1
	}
	p.value<-k/nrepl
}
else
{
	p.value<-2*(1-pnorm(abs(t/(sqrt(((gamma(1+p))^(2/p))*(-1-(1/(p^2))+(gamma(1+2*p))/((p^2)*(gamma(1+p))^2)))))))
}
RVAL<-list(statistic=c(T=t), p.value=p.value, method="Atkinson test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}