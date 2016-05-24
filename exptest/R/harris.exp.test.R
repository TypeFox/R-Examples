harris.exp.test<-function(x, R=length(x)/4, simulate.p.value=FALSE, nrepl=2000)
{
DNAME <- deparse(substitute(x))
R=round(R)
n<-length(x)
x<-sort(x)
x<-c(0,x)
D<-(n:1)*(x[2:(n+1)]-x[1:n])
t<-((sum(D[1:R])+sum(D[(n-R+1):n]))/(2*R))/((sum(D[(R+1):(n-R)]))/(n-2*R))
l<-0
if(simulate.p.value)
{
	for(i in 1:nrepl)
	{
	z<-rexp(n)
	z<-sort(z)
	z<-c(0,z)
	D<-(n:1)*(z[2:(n+1)]-z[1:n])
	T<-((sum(D[1:R])+sum(D[(n-R+1):n]))/(2*R))/((sum(D[(R+1):(n-R)]))/(n-2*R))
	if(abs(T)>abs(t)) l=l+1
	}
	p.value<-l/nrepl
}
else
{	
	p.value<-2*min(pf(t,4*R,2*(n-2*R)),1-pf(t,4*R,2*(n-2*R)))
}
RVAL<-list(statistic=c(Q=t), p.value=p.value, method="Harris modification of Gnedenko's F-test",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}