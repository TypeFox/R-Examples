neyman.unif.test<-function(x, nrepl=2000,k=5)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
po<-legendre.polynomials(k,normalized=TRUE)
ne<-0
for(i in 1:k)
{
	v<-(sum(sqrt(2)*polynomial.values(po[i+1],2*x-1)[[1]])/sqrt(n))^2
	ne<-ne+v
}

for(i in 1:nrepl)
{
	z<-runif(n)
	z<-sort(z)
	Ne<-0
	for(j in 1:k)
	{
		V<-(sum(sqrt(2)*polynomial.values(po[j+1],2*z-1)[[1]])/sqrt(n))^2
		Ne<-Ne+V
	}
	if (Ne>ne) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(N=ne), p.value=p.value, method="Neyman test for uniformity",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}