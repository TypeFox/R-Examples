kimber.exp.test<-function(x, nrepl=2000)
{
DNAME <- deparse(substitute(x))
l<-0
n<-length(x)
x<-sort(x)
y<-mean(x)
s<-(2/pi)*asin(sqrt(1-exp(-(x/y))))
r<-(2/pi)*asin(sqrt((c(1:n) - 0.5)/n))
d<-max(abs(r-s))
for(i in 1:nrepl)
{
z<-rexp(n)
z<-sort(z)
y<-mean(z)
s<-(2/pi)*asin(sqrt(1-exp(-(z/y))))
r<-(2/pi)*asin(sqrt((c(1:n) - 0.5)/n))
D<-max(abs(r-s))
if (D>d) l=l+1
}
p.value<-l/nrepl
RVAL<-list(statistic=c(D=d), p.value=p.value, method="Kimber-Michael test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}