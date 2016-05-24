# expected value
Edweibull3<-Vectorize(function(c, beta, eps=0.0001)
{
if(beta==0)
{E<-exp(-c)/(1-exp(-c))}
else
{
nmax<-qdweibull3(1-eps, c, beta)*2
E<-0
for(i in 1:nmax)
{
E<-E+exp(-c*(sum((1:i)^beta)))
}
}
E
}
)
# second order moment
E2dweibull3<-Vectorize(function(c, beta, eps=0.0001)
{
if(beta==0)
{E<-exp(-c)*(1+exp(-c))/(1-exp(-c))^2}
else
{
nmax<-qdweibull3(1-eps, c, beta)*2
E<-0
for(i in 1:nmax)
{
E<-E+2*i*exp(-c*(sum((1:i)^beta)))
}
E-Edweibull3(c,beta)
}
}
)