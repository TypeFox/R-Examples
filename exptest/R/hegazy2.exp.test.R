hegazy2.exp.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
b<--log(1-c(1:n)/(n+1))
t<-(n^(-1))*sum((x-b)^2)
for(i in 1:nrepl)
{
z<-rexp(n)
z<-sort(z)
T<-(n^(-1))*sum((z-b)^2)
if (T>t) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(T=t), p.value=p.value, method="Hegazy-Green test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}