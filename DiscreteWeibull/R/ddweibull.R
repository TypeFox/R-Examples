ddweibull<-Vectorize(function(x, q, beta, zero=FALSE)
{
if(x!=floor(x))
{stop("Non-integer x = ",x)}
if(zero)
q^x^beta-q^(x+1)^beta
else
q^(x-1)^beta-q^x^beta
}
)

pdweibull<-Vectorize(function(x, q, beta, zero=FALSE)
{
if (zero)
1-q^(x+1)^beta
else
1-q^x^beta
}
)

qdweibull<-Vectorize(function(p, q, beta, zero=FALSE)
{
if(zero)
ceiling((log(1-p)/log(q))^(1/beta)-1)
else
ceiling((log(1-p)/log(q))^(1/beta))
}
)

rdweibull<-function(n, q, beta, zero=FALSE)
{
u<-runif(n)
if(zero)
{s<-ceiling((log(1-u)/log(q))^(1/beta))-1}
else
{s<-ceiling((log(1-u)/log(q))^(1/beta))}
s
}