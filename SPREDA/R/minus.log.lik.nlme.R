minus.log.lik.nlme <-
function(dat, coef, random.eff)
{
  dfs=dat$dfs
  lam=dat$lam
  dat=dat$dat
  A=coef[1]
  B=coef[2]
  beta=coef[3:(length(coef)-4)]
  ss0=coef[length(coef)-3]
  ss1=coef[length(coef)-2]
  rho=coef[length(coef)-1]
  ss=coef[length(coef)]
  Sigma.D=matrix(c(ss0^2,rho*ss0*ss1,rho*ss0*ss1,ss1^2),2,2)/(ss^2)
  ids=unique(dat[,1])
  nn=length(ids)
  res=0
  
  for(i in 1:nn){
    idx=(dat[,1]==ids[i])
    n.idx=sum(idx)
    dat.i=dat[idx,]
    y=dat.i[, 3]
    Dmat=as.matrix(dat.i[, 4:(length(beta)+3)])
    rand=random.eff[i,]
    y.hat=-as.matrix(A*exp(rand[1])/(1+exp(-log(Dmat%*%beta)/(B*exp(rand[2])))))
    Zhat=matrix(0, nrow=n.idx, ncol=2)
    Zhat[,1]=y.hat
    Zhat[,2]=y.hat*(1+exp(-log(Dmat%*%beta)/(B*exp(rand[2]))))^(-1)*(-1)*exp(-log(Dmat%*%beta)/(B*exp(rand[2])))*(log(Dmat%*%beta)/(B*exp(rand[2])))
    Vhat=diag(n.idx)+Zhat%*%Sigma.D%*%t(Zhat)
    resi=as.matrix(y-y.hat+Zhat%*%rand)
    tmp=-0.5*log(det(ss^2*Vhat))-0.5*ss^(-2)*(t(resi)%*%solve(Vhat)%*%resi)
    res=res+tmp
    
  }
  return(-res)
  
}
