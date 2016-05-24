ks.exp.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
y<-x/mean(x)
z<-sort(1-exp(-y))
j1<-c(1:n)/n
m1<-max(j1-z)
j2<-(c(1:n)-1)/n
m2<-max(z-j2)
t<-max(m1,m2)
for(i in 1:nrepl)
{
s<-rexp(n)
y<-s/mean(s)
z<-sort(1-exp(-y))
m1<-max(j1-z)
m2<-max(z-j2)
T<-max(m1,m2)
if (T>t) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(KSn=t), p.value=p.value, method="Kolmogorov-Smirnov test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}