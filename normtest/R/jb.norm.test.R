jb.norm.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x1<-sum(x)/n
b1<-sqrt(n)*(sum((x-x1)^3))/((sum((x-x1)^2))^(3/2))
b2<-n*(sum((x-x1)^4))/((sum((x-x1)^2))^2)
t<-(n/6)*((b1)^2 + ((b2-3)^2)/4)
for(i in 1:nrepl)
{
	z<-rnorm(n)
	z1<-sum(z)/n
	a1<-sqrt(n)*(sum((z-z1)^3))/((sum((z-z1)^2))^(3/2))	
	a2<-n*(sum((z-z1)^4))/((sum((z-z1)^2))^2)
	T<-(n/6)*((a1)^2 + ((a2-3)^2)/4)
	if (T>t) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(JB=t), p.value=p.value, method="Jarque-Bera test for normality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}