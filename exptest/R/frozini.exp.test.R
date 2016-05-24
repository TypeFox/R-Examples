frozini.exp.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
y<-mean(x)
b<-(1/sqrt(n))*sum(abs(1-exp(-x/y)-(c(1:n)-0.5)/n))
for(i in 1:nrepl)
{
	z<-rexp(n)
	z<-sort(z)
	y<-mean(z)
	B<-(1/sqrt(n))*sum(abs(1-exp(-z/y)-(c(1:n)-0.5)/n))
	if (B>b) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(B=b), p.value=p.value, method="Frozini test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}