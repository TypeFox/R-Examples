pietra.exp.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
xm<-mean(x)
t<-sum(abs(x-xm))/(2*n*xm)
for(i in 1:nrepl)
{
s<-rexp(n)
sm<-mean(s)
T<-sum(abs(s-sm))/(2*n*sm)
l=l+(T>t)
}
p.value<-2*min(l/nrepl,1-l/nrepl)
RVAL<-list(statistic=c(Pn=t), p.value=p.value, method="Test for exponentiality based on the Pietra statistic",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}