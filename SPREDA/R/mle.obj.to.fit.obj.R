mle.obj.to.fit.obj <-
function(obj)
{
  random.names=obj$dat$random.names
  dat=obj$dat$dat
  pars=obj$coef
  n.random=length(random.names)
  n.dat=ncol(dat)
  n.pars=length(pars)
  random.index=match(random.names, colnames(dat), 0L)
  if(n.pars==4){
    ss0=exp(pars[1])
    ss1=exp(pars[2])
    rho=2/(1+exp(-pars[3]))-1
    ss=exp(pars[4])
    Sigma=matrix(c(ss0^2,rho*ss0*ss1,rho*ss0*ss1,ss1^2),2,2)
    ss.vec=c(ss0, ss1, rho, ss)
  }else if(n.pars==2){
    ss.r=exp(pars[1])
    ss=exp(pars[2])
    Sigma=ss.r^2
    ss.vec=c(ss.r, ss)
  }
  
  
  ids=unique(dat[,1])
  mm=length(ids)
  yy=dat[,n.dat]
  error=yy
  fitted.val=error
  ran.eff=matrix(0,nrow=mm,ncol=n.random)
  
  
  
  for(j in 1:mm)
  {
    id.x=(dat[,1]==ids[j])
    tmp1=yy[id.x]
    Z=dat[id.x,random.index]
    Z=as.matrix(Z)
    ww=dim(Z)[1]
    SS=Z%*%Sigma%*%t(Z)+ss^2*diag(ww)
    SS.inv=solve(SS)
    ran.val=Sigma%*%t(Z)%*%SS.inv%*%tmp1
    ran.eff[j,]=ran.val
    t.val=Z%*%ran.val
    fitted.val[id.x]=t.val
    error[id.x]=tmp1-t.val
  }
  
  random=list(ran.eff)
  coef=list(random=random)
  
  res=list(ss.vec=ss.vec, fitted=fitted.val,error=error,coef=coef)
  
  return(res)
}
