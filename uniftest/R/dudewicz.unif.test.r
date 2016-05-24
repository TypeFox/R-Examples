dudewicz.unif.test<-function(x, nrepl=2000,m=length(x)/2)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
a<-x[1]
b<-x[n]
for(i in 1:m) x<-c(a,x,b)
h<-(-1)/n*sum(log(n/2/m*(x[c(1:n)+2*m]-x[c(1:n)]),2))
for(i in 1:nrepl)
{
	z<-runif(n)
	z<-sort(z)
	a<-z[1]
	b<-z[n]
	for(i in 1:m) z<-c(a,z,b)
	H<-(-1)/n*sum(log(n/2/m*(z[c(1:n)+2*m]-z[c(1:n)]),2))
	if (H>h) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(H=h), p.value=p.value, method="Dudewicz-van der Meulen test for uniformity",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}