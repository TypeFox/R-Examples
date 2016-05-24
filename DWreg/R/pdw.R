pdw<-function(x,q=exp(-1),beta=1)
{
if(x<0)
res<-0
else
res<-1-q^((floor(x)+1)^(beta))
return(res)
}