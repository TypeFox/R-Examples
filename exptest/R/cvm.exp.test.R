cvm.exp.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
y<-x/mean(x)
z<-sort(1-exp(-y))
c<-(2*c(1:n)-1)/(2*n)
z<-(z-c)^2
t<-1/(12*n)+sum(z)
for(i in 1:nrepl)
{
s<-rexp(n)
y<-mean(s)
z<-sort(1-exp(-y))
c<-(2*c(1:n)-1)/(2*n)
z<-(z-c)^2
T<-1/(12*n)+sum(z)
if (T>t) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(Wn=t), p.value=p.value, method="Cramer-von Mises test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}