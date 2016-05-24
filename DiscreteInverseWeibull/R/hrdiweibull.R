hrdiweibull<-Vectorize(function(x, q, beta)
{
if(x!=floor(x) | x<1)
{stop("Non-feasible value of x!",x)}
(q^x^(-beta)-q^(x-1)^(-beta))/(1-q^(x-1)^(-beta))
}
)

