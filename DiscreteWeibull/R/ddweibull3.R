# random sample
rdweibull3<-function(n, c, beta)
{
u<-runif(n)
qdweibull3(u, c, beta)
}

# probability mass function
ddweibull3<-Vectorize(function(x, c, beta)
{
if(x==0)
d<-1-exp(-c)
else
{
s<-1:x
d<-(1-exp(-c*(x+1)^beta))*exp(-c*(sum(s^beta)))
}
d
}
)
# cumulative distribution function
pdweibull3<-Vectorize(function(x, c, beta)
{
if(x==0)
1-exp(-c)
else
{
s<-1:(x+1)
1-exp(-c*sum(s^beta))
}
}
)
# inverse of the cumulative distribution function
qdweibull3<-Vectorize(function(p, c, beta)
{
if(p<=1-exp(-c))
q<-0
else
{
q<-1
while(p>1-exp(-c*sum((1:(q+1))^beta)))
q<-q+1
}
q
}
)
######################################
hdweibull3<-function(x,c,beta)
{
1-exp(-c(x+1)^beta)
}


