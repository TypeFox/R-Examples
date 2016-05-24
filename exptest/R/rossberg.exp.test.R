rossberg.exp.test<-function(x)
{
DNAME <- deparse(substitute(x))
n<-length(x)
s<-0
sh<-0
sg<-0
for (m in 1:n)
{
h<-0
for (i in 1:(n-2))
for (j in (i+1):(n-1))
for (k in (j+1):n)
{
if ((x[i]+x[j]+x[k]-2*min(x[i],x[j],x[k])-max(x[i],x[j],x[k]))<x[m])
{
h=h+1
}
}
h=((6*factorial(n-3))/factorial(n))*h
sh=sh+h
}
for (m in 1:n)
{
g<-0
for (i in 1:(n-1))
for (j in (i+1):n)
{
if (min(x[i],x[j])<x[m])
{
g=g+1
}
}
g=((2*factorial(n-2))/factorial(n))*g
sg=sg+g
}
s=sh-sg
s<-s/n
v<-sqrt(n)*abs(s)/sqrt(52/1125)
p.value<-2*(1-pnorm(v))
RVAL<-list(statistic=c(Sn=s), p.value=p.value, method="Test for exponentiality based on Rossberg's characterization",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}