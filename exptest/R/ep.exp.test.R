ep.exp.test<-function(x, simulate.p.value=FALSE, nrepl=2000)
{
DNAME <- deparse(substitute(x))
n<-length(x)
y<-x/mean(x)
ep<-sqrt(48*n)*sum(exp(-y)-1/2)/n
l<-0
if(simulate.p.value)
{
	for(i in 1:nrepl)
	{
	z<-rexp(n)
	y<-z/mean(z)
	EP<-sqrt(48*n)*sum(exp(-y)-1/2)/n
	if (abs(EP)>abs(ep)) l=l+1
	}
	p.value<-l/nrepl
}
else
{
	p.value<-2*(1-pnorm(abs(ep)))
}
RVAL<-list(statistic=c(EPn=ep), p.value=p.value, method="The test for exponentiality of Epps and Pulley",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}