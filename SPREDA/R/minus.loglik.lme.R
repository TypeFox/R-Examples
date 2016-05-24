minus.loglik.lme <-
function(dat, pars)
{
  random.names=dat$random.names
  newdat=dat$dat
  
  n.dat=ncol(newdat)
  n.pars=length(pars)
  if(n.pars==4){
    ss0=exp(pars[1])
    ss1=exp(pars[2])
    rho=2/(1+exp(-pars[3]))-1
    ss=exp(pars[4])
    Sigma=matrix(c(ss0^2,rho*ss0*ss1,rho*ss0*ss1,ss1^2),2,2)
    random.index=match(random.names, colnames(newdat), 0L)
  }else if(n.pars==2){
    ss.r=exp(pars[1])
    ss=exp(pars)[2]
    Sigma=ss.r^2
    random.index=match(random.names, colnames(newdat), 0L)
  }
  
  
  ids=unique(newdat[,1])
  mm=length(ids)
  error=newdat[,n.dat]
  
  loglik=0
  
  for(j in 1:mm)
  {
    id.x=(newdat[,1]==ids[j])
    tmp1=error[id.x]
    Z=newdat[id.x, random.index]
    Z=as.matrix(Z)
    ww=dim(Z)[1]
    SS=Z%*%Sigma%*%t(Z)+ss^2*diag(ww)
    SS.inv=solve(SS)
    ll=-(ww/2)*log(2*pi)-.5*log(det(SS))-.5*as.vector(t(tmp1)%*%SS.inv%*%tmp1)
    loglik=loglik+ll
  }
  
  res=(-1)*loglik

  return(res)
  
}
