lorenz.exp.test<-function(x, p=0.5, simulate.p.value=FALSE, nrepl=2000)
{
DNAME <- deparse(substitute(x))
k<-0
n<-length(x)
x<-sort(x)
l<-sum(x[1:floor(n*p)])/(sum(x))
if (simulate.p.value)
{
	for(i in 1:nrepl)
	{
	z<-rexp(n)
	z<-sort(z)
	L<-sum(z[1:floor(n*p)])/(sum(z))
	if (abs(L)>abs(l)) k=k+1
	}
	p.value<-k/nrepl
}
else
{
	p.value<-2*(1-pnorm(abs(sqrt(n)*(l-p-(1-p)*log(1-p)))))
}
RVAL<-list(statistic=c(L=l), p.value=p.value, method="Lorenz test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}