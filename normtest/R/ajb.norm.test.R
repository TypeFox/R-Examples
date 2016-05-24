ajb.norm.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x1<-sum(x)/n
b1<-sqrt(n)*(sum((x-x1)^3))/((sum((x-x1)^2))^(3/2))
b2<-n*(sum((x-x1)^4))/((sum((x-x1)^2))^(2))
Eb2<-3*(n-1)/(n+1)
varb2<-(24*n*(n-2)*(n-3))/((n+1)*(n+1)*(n+3)*(n+5))
varb1<-(6*(n-2))/((n+1)*(n+3))
t<-(b1^2)/varb1 + (b2-Eb2)^2/varb2
for(i in 1:nrepl)
{
	z<-rnorm(n)
	z1<-sum(z)/n
	a1<-sqrt(n)*(sum((z-z1)^3))/((sum((z-z1)^2))^(3/2))	
	a2<-n*(sum((z-z1)^4))/((sum((z-z1)^2))^(2))
	T<-(a1^2)/varb1 + (a2-Eb2)^2/varb2
	if (abs(T)>abs(t)) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(AJB=t), p.value=p.value, method="Adjusted Jarque-Bera test for normality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}