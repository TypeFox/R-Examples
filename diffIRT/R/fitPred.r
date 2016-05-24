fitPred=function(ai,vi,sd2A,sd2V,order,model=1,A,W){     # calculates predicted values up to a certain order for D diffusion and Q diffusion
  nit=length(ai)
  nq=length(A)
  X=fitItems(nit,order)
  N=fitCumChoose(nit,order)
  ou=matrix(N,nit)
  vi=matrix(vi)
  ai=matrix(ai)
  V=A
  I=matrix(1,nq)
  V=matrix(V*sqrt(sd2V))%x%I
  A=I%x%matrix(A*sqrt(sd2A))
  W=(matrix(W)%x%I) * (I%x% matrix(W))

  if(model==1) dev=kronecker(V,t(vi),"-")                   #D diffusion
  if(model==2) dev=kronecker(exp(V),exp(t(vi)),"/")         #Q diffusion

  amat=kronecker(exp(A),exp(t(ai)),"/")
  dev2=amat*dev
  p=1/(1+exp(-dev2))

  ou=matrix(,nrow(X),1)
  for(i in 1:nrow(X)) ou[i,1]=sum(W*apply(t(t(p)^(X[i,]))*t(t(1-p)^(1-X[i,])),1,prod))
  return(ou)
}
