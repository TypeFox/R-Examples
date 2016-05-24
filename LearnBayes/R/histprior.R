histprior=function(p,midpts,prob)
{
binwidth=midpts[2]-midpts[1]
lo=round(10000*(midpts-binwidth/2))/10000
val=0*p
for (i in 1:length(p))
{
 val[i]=prob[sum(p[i]>=lo)]
}
return(val)
}