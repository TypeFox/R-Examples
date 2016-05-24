miniusloglik.sev <-
function(dat, pars)
{
 #column 1 event time
 #column 2 for event indicator
 tt=dat[,1]
 delta=dat[,2]
 mu=pars[1]
 sigma=exp(pars[2])
 zz=(log(tt)-mu)/sigma
 ff=dsev(zz)/(sigma*tt)
 FF=psev(zz)
 
 ll=delta*log(ff)+(1-delta)*log(1-FF)
 res=(-1)*sum(ll)
 return(res)
}
