Ediweibull <-
function(q, beta, eps=1e-04, nmax=1000)
{
xmax<-max(2*qdiweibull(1-eps, q, beta), nmax)

x<-1:xmax
if(beta>1)
{
EX=1+sum(1-q^x^(-beta))
}
else
{EX=Inf}
if(beta>2)
{
EX2=2*sum(x*(1-q^x^(-beta)))+EX
}
else
{EX2=Inf}
list(EX=EX,EX2=EX2)
}

