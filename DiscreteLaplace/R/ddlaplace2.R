ddlaplace2 <-
function(x,p,q)
{
palaplace2(x+1,p,q)-palaplace2(x,p,q)
}


palaplace2 <-Vectorize(
function(x,p,q)
{
if(x>0)
1-log(q)/log(p*q)*p^x
else if(x<0)
log(p)/log(p*q)*q^(-x)
else if (x==0)
log(p)/log(p*q)
}
)

pdlaplace2 <-
function(x,p,q)
{
palaplace2(x+1,p,q)
}


qdlaplace2 <-Vectorize(
function(prob,p,q)
{
h<-1-log(q)/log(p*q)*p
if(prob>=h)
{
ceiling((log(1-prob)+log(log(p*q)/log(q)))/log(p))-1
}
else
{
ceiling((log(log(p)/log(p*q))-log(prob))/log(q))-1
}
}
)

rdlaplace2 <-
function(n,p,q)
{
u<-runif(n)
qdlaplace2(u,p,q)
}
