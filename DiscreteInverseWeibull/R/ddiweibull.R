ddiweibull<-Vectorize(function(x, q, beta)
{
if(x!=floor(x) | x<=0)
{stop("Non-feasible value for x = ",x)}
if(x>1)
{
q^(x)^(-beta)-q^(x-1)^(-beta)
}
else
{
q
}
}
)

pdiweibull<-Vectorize(function(x, q, beta)
{
if(x<1)
{
0
}
else
{
q^floor(x)^(-beta)
}
}
)

qdiweibull<-Vectorize(function(p, q, beta)
{
if(p>=1 | p<=0)
{
stop("Please choose 0<p<1!")
}
else
{
ceiling((log(p)/log(q))^(-1/beta))
}
}
)

rdiweibull <-
function(n, q, beta)
{
u<-runif(n)
ceiling((log(u)/log(q))^(-1/beta))
}


