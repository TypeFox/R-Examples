ww.exp.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
t<-max(x)/min(x)
for(i in 1:nrepl)
{
s<-rexp(n)
T<-max(s)/min(s)
if (T>t) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(Q=t), p.value=p.value, method="Wong and Wong's test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}