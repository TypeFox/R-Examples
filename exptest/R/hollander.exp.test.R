hollander.exp.test<-function(x, simulate.p.value=FALSE, nrepl=2000)
{
DNAME <- deparse(substitute(x))
n<-length(x)
t<-0
for (i in 1:n)
for (j in 1:n)
for (k in 1:n)
{
	if ((i!=j) & (i!=k) & (j<k) & (x[i]>x[j]+x[k]))
	{
		t=t+1
	}
}
t<-(2/(n*(n-1)*(n-2)))*t
if(simulate.p.value)
{
	l<-0
	for (m in 1:nrepl)
	{	
		z<-rexp(n)
		T<-0
		for (i in 1:n)
		for (j in 1:n)
		for (k in 1:n)
		{
			if ((i!=j) & (i!=k) & (j<k) & (z[i]>z[j]+z[k]))
			{
				T=T+1
			}
		}			
		T<-(2/(n*(n-1)*(n-2)))*T
		if (abs(T)>abs(t))
		{
			l<-l+1
		}
	}
	p.value<-l/nrepl
}
else
{
	p.value<-2*(1-pnorm(sqrt(n)*abs(t-1/4)/sqrt(5/432)))
}
RVAL<-list(statistic=c(T=t), p.value=p.value, method="Hollander-Proshan test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}