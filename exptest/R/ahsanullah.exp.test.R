ahsanullah.exp.test<-function(x, simulate.p.value=FALSE, nrepl=2000)
{
DNAME <- deparse(substitute(x))
n<-length(x)
h<-0
g<-0
for (k in 1:n)
{
	for (i in 1:n)
	for (j in 1:n)
	{
		if (abs(x[i]-x[j])<x[k])
		{
			h=h+1
		}
		if (2*min(x[i],x[j])<x[k])
		{
			g=g+1
		}
	}
}
r<-(h-g)/(n^3)
v<-sqrt(n)*abs(r)/sqrt(647/4725)
if(simulate.p.value)
{
	l<-0
	for(m in 1:nrepl)
	{
	z<-rexp(n)
	h<-0
	g<-0
	for (k in 1:n)
	{
		for (i in 1:n)
		for (j in 1:n)
		{
			if (abs(z[i]-z[j])<z[k])
			{
				h=h+1
			}
			if (2*min(z[i],z[j])<z[k])
			{
				g=g+1
			}
		}
	}
	R<-(h-g)/(n^3)
	if (abs(R)>abs(r))
	{
		l<-l+1
	}
	}
	p.value<-l/nrepl
}
else
{
	p.value<-2*(1-pnorm(v))
}
RVAL<-list(statistic=c(In=r), p.value=p.value, method="Test for exponentiality based on Ahsanullah's characterization",data.name = DNAME)
class(RVAL)<-"htest"
return(RVAL)
}