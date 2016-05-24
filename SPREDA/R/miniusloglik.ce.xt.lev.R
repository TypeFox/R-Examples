miniusloglik.ce.xt.lev <-
function(dat, pars)
{
  mu=pars[1] #mu
  npars=length(pars)
  sigma=exp(pars[npars])   #sigma
  beta.vec=pars[2:(npars-1)] #beta: coefficient for covariates
  
  failure.dat=dat$failure.dat
  aux.inf=dat$aux.inf
  xt.obj=dat$xt.obj
  
  wts.mat=aux.inf$wts.mat
  npts.vec=aux.inf$npts.vec
  
  x.val=xt.obj$x.val
  npts=xt.obj$npts
  n=xt.obj$n
  
  para.vec=kronecker(beta.vec,rep(1,npts))
  para.mat=kronecker(t(para.vec),rep(1,n))
  
  beta.xt.mat=para.mat*x.val
  tmp.mat=kronecker(rep(1,npars-2),diag(npts))
  
  sum.beta.xt.mat=beta.xt.mat%*%tmp.mat
  sum.beta.xt.mat=exp(sum.beta.xt.mat)
  
  int.g.beta.xt=rowSums(wts.mat*sum.beta.xt.mat)
  
  g.beta.xt=as.vector(sum.beta.xt.mat)[1:n+n*(npts.vec-1)]
  
  delta=failure.dat[,"delta"]
  zz=(log(int.g.beta.xt)-mu)/sigma
  ff=dlev(zz)/(int.g.beta.xt*sigma)
  FF=plev(zz)
  
  ll=delta*log(g.beta.xt*ff)+(1-delta)*log(1-FF)
  res=(-1)*sum(ll)
  
  return(res)
}
