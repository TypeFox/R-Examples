we.exp.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
m<-mean(x)
v<-var(x)
t<-(n-1)*v/(n^2*m^2)
for(i in 1:nrepl)
{
s<-rexp(n)
m<-mean(s)
v<-var(s)
T<-(n-1)*v/(n^2*m^2)
if (T>t) l=l+1
}
p.value<-2*min(l/nrepl,1-l/nrepl)
RVAL<-list(statistic=c(WE=t), p.value=p.value, method="WE test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}