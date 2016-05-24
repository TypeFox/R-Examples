shapiro.exp.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
y<-mean(x)
w<-n*(y-x[1])^2/((n-1)*sum((x-y)^2))
for(i in 1:nrepl)
{
s<-rexp(n)
s<-sort(s)
y<-mean(s)
W<-n*(y-s[1])^2/((n-1)*sum((s-y)^2))
if (W<w) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(W=w), p.value=p.value, method="Shapiro-Wilk test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}