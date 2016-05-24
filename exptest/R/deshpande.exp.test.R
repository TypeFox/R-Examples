deshpande.exp.test<-function(x, b=0.44, simulate.p.value=FALSE, nrepl=2000)
{
DNAME <- deparse(substitute(x))
n<-length(x)
j<-0
v<-0
	for (i in 1:n)
	for (k in 1:n)
	{
		if ((i!=k) & (x[i]>b*x[k]))
		{
			j=j+1
		}
	}
j<-j/(n*(n-1))
v<-1 + b/(b+2) + 1/(2*b+1) + (2*(1-b))/(b+1) - (2*b)/(b*b + b +1) - 4/((b+1)*(b+1))
if(simulate.p.value)
{
	l<-0
	for (m in 1:nrepl)
	{
		z<-rexp(n)
		J<-0
		for (i in 1:n)
		for (k in 1:n)
		{
			if ((i!=k) & (z[i]>b*z[k]))
			{
				J=J+1
			}
		}
		J<-J/(n*(n-1))
		if (abs(J)>abs(j))
		{
			l<-l+1
		}
	}
	p.value<-l/nrepl
}
else
{
	p.value<-2*(1-pnorm(sqrt(n)*abs(j-1/(b+1))/sqrt(v)))
}
RVAL<-list(statistic=c(J=j), p.value=p.value, method="Deshpande test for exponentiality",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}