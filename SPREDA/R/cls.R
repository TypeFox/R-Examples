cls <-
function(y,X)
{
  y0=y
  X0=X
  
  y=Px(X)%*%y
  nn=dim(X)[1]
  p=dim(X)[2]
  Z=matrix(0,nrow=nn,ncol=p)
  for(i in 1:p)
  {
    Z[,i]=-(diag(nn)-Px(X[,-i]))%*%X[,i]
  }
  
  ind.X=rep(1,p)
  ind.Z=rep(2,p)
  
  bb=solve(t(X)%*%X)%*%t(X)%*%y
  bb=as.vector(bb)
  bb0=rep(max(abs(bb)),p)
  zz=as.vector(X%*%bb0)
  
  #print(bb)
  while(any(bb<0))
  {
    #print(bb)
    tt=bb0/(bb0-bb)
    id=min((1:p)[tt==min(tt[bb<0])])
    tmp=X[,id]
    X[,id]=Z[,id]
    Z[,id]=tmp
    
    tmp=ind.X[id]
    ind.X[id]=ind.Z[id]
    ind.Z[id]=tmp
    
    zz=zz+min(tt[id],1)*(y-zz)
    bb0=as.vector(solve(t(X)%*%X)%*%t(X)%*%zz)
    bb=as.vector(solve(t(X)%*%X)%*%t(X)%*%y)
    #print(bb)
    # print(ind.X)
  }
  bb[ind.X==2]=0
  
  betahat=bb
  yhat=as.vector(X%*%bb)
  
  #print(betahat)
  res=list(y=y0,X=X0,betahat=betahat,yhat=yhat)
  return(res)
}
