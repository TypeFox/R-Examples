miniusloglik.logis.wts <-
function(dat, pars)
{
  n=nrow(dat)
  index=match("time", colnames(dat), 0)
  if(index==1){
    dat=as.data.frame(cbind(start=rep(0, n), stop=dat[,1], status=dat[,2], dat[,-c(1,2)]))
  }
  tt=dat[,"stop"]
  delta=dat[,"status"]
  wts=dat[,"wts"]
  cov=as.matrix(dat[,-c(1:4)])
  
  m=length(pars)
  #mu=pars[1]
  sigma=exp(pars[m])
  beta=as.matrix(pars[-m])
  zz=(log(tt)-cov%*%beta)/sigma
  ff=dlogis(zz)/(sigma*tt)
  FF=plogis(zz)
  
  tt0=dat[,"start"]
  FF0=rep(0, n)
  idx=(tt0==0)
  if(sum(idx)>0){
    tt0.trun=tt0[!idx]
    zz0=(log(tt0.trun)-cov[!idx,]%*%beta)/sigma
    FF0[!idx]=plogis(zz0)
  }
  
  ll=delta*log(ff/(1-FF0))+(1-delta)*log((1-FF)/(1-FF0))
  res=(-1)*sum(wts*ll)
  return(res)
}
