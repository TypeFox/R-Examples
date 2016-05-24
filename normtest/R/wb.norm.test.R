wb.norm.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
m<-qnorm((c(1:n)-(3/8))/(n+(1/4)))
c<-m/sqrt(sum(m*m))
x<-sort(x)
s<-var(x)*(n-1)
wb<-((sum(c*x))^2)/s
for(i in 1:nrepl)
{
	z<-rnorm(n)
	z<-sort(z)
	s<-var(z)*(n-1)
	Wb<-((sum(c*z))^2)/s
	if (Wb<wb) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(WB=wb), p.value=p.value, method="Weisberg-Bingham test for normality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}