quesenberry.unif.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-c(0,sort(x),1)
q<-sum((x[c(2:(n+2))]-x[c(1:(n+1))])^2)+sum((x[c(2:(n+1))]-x[c(1:n)])*(x[c(3:(n+2))]-x[c(2:(n+1))]))
for(i in 1:nrepl)
{
	z<-runif(n)
	z<-c(0,sort(z),1)
	Q<-sum((z[c(2:(n+2))]-z[c(1:(n+1))])^2)+sum((z[c(2:(n+1))]-z[c(1:n)])*(z[c(3:(n+2))]-z[c(2:(n+1))]))
	if (Q>q) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(Q=q), p.value=p.value, method="Greenwood-Quesenberry-Miller test for uniformity",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}