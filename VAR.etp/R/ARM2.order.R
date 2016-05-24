ARM2.order <-
function(x,y,pmax)
{
aicmat=matrix(NA,nrow=pmax);
bicmat=matrix(NA,nrow=pmax);

for (i in 1:pmax)
{
M = ARM2(x,y,p=i)
aicmat[i,]=M$aic;bicmat[i,]=M$bic;
}
p.aic=which.min(aicmat)
p.bic=which.min(bicmat)
return(list(p.aic=p.aic,p.bic=p.bic))
}
