Confint.mean.two_stage <- function(L, alpha, presample, addsample)
{
n0=length(presample)
s0=sd(presample)
z=(L/2)^2/qt(1-alpha/2,n0-1)^2
n=n0+length(addsample)
a=z*(s0^2)
c=1/(n-n0)-a
b=-2/(n-n0)
d=(n/n0)/(n-n0)
if (b^2 > 4*c*d)
{a01=n0/n-sqrt((b^2-4*c*d)/(4*d^2))
 a02=n0/n+sqrt((b^2-4*c*d)/(4*d^2))}
else {a01=n0/n
      a02=(n-n0)/n}
a11=(1-a01)
a12=(1-a02)
ln1=a01*mean(presample)+a11*mean(addsample)
ln2=a02*mean(presample)+a12*mean(addsample)
if ((ln1-L/2<mean(presample)) & (ln1+L/2>mean(presample)))
    {ln=ln1}
else {ln=ln2}
list("Significance level"=alpha,
"Length of confidence interval"=L,
"Length of presample"=n0,
"Number of additional observations"=n-n0,
"Total number of observations"=n,
"confidence interval"=c(ln-L/2,ln+L/2))
}

#pres=rnorm(16,10,1.2)
#v=add_size.mean.two_stage(0.9,0.05,pres)$"Number of additional observations"
#confint.mean.two_stage(0.9,0.05,pres,rnorm(v,10,1.2))
