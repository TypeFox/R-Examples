kochar.exp.test<-function(x, simulate.p.value=FALSE, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
u<-c(1:n)/(n+1)
j<-2*(1-u)*(1-log(1-u)) - 1
t<-sqrt(108*n/17)*(sum(j*x))/(sum(x))
if (simulate.p.value)
{
	for(i in 1:nrepl)
	{
	z<-rexp(n)
	z<-sort(z)
	T<-sqrt(108*n/17)*(sum(j*z))/(sum(z))
	if (abs(T)>abs(t)) l=l+1
	}
	p.value<-l/nrepl
}
else
{
	p.value<-2*(1-pnorm(abs(t)))
}
RVAL<-list(statistic=c(T=t), p.value=p.value, method="Kochar test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}