fitDeriv=function(ai,vi,sd2A,sd2V,order,model=1,A,W){
  nq=length(A)
  nit=length(ai)
  na=fitCumChoose(nit,order)+1
  ret=matrix(,na,2*nit+2)
  X=fitItems(nit,order)
  ai=matrix(exp(ai))
  vi=matrix(vi)
  if(toupper(model)=="Q") vi=exp(vi)
  V=A
  I=matrix(1,nq)
  V=matrix(V*sqrt(sd2V))%x%I
  A=I%x%matrix(A*sqrt(sd2A))
  W=(matrix(W)%x%I) * (I%x% matrix(W))

  amat=kronecker(exp(A),t(ai),"/")                # necessary for both Q and D diffusion
  amat_sq=kronecker(exp(A),t(ai)^2,"/")

  if(model==1) dev=kronecker(V,t(vi),"-")         #D diffusion
  if(model==2){
                       dev=kronecker(exp(V),t(vi),"/")         #Q diffusion
                       dev_sq=kronecker(exp(V),t(vi)^2,"/")
                       logit_v_sq=amat*dev_sq
                       }
logit=amat*dev
logit_a_sq=amat_sq*dev
logit[logit<(-709)]=-709   #prevent overflow

p=1/(1+exp(-logit))
if(model==1){
  derPai=-p^2*exp(-logit)*logit_a_sq
  derPvi=-p^2*exp(-logit)*amat
  }
if(model==2){
  derPai=-p^2*exp(-logit)*logit_a_sq
  derPvi=-p^2*exp(-logit)*logit_v_sq
  }
X2=2*X-1

  ou=matrix(,nrow(X),nit)
if(nq>1){
  for(j in 1:na){
    ret[j,2*nit+1]=sum(W*apply((t(t(p)^X[j,])*t(t(1-p)^(1-X[j,]))),1,prod)*(-1/sqrt(sd2A)*(1-A^2/sd2A)))  #derivatives to sd(A) and sd(V)
    ret[j,2*nit+2]=sum(W*apply((t(t(p)^X[j,])*t(t(1-p)^(1-X[j,]))),1,prod)*(-1/sqrt(sd2V)*(1-V^2/sd2V)))
  for(i in 1:nit){
    ret[j,i]= sum(W*derPai[,i]*X2[j,i]*apply((t(t(p)^X[j,])*t(t(1-p)^(1-X[j,])))[,-i],1,prod))
    ret[j,i+nit]= sum(W*derPvi[,i]*X2[j,i]*apply((t(t(p)^X[j,])*t(t(1-p)^(1-X[j,])))[,-i],1,prod))
  }
  }
}
if(nq==1){
  for(j in 1:na){
    ret[j,2*nit+1]=sum(W*prod((t(t(p)^X[j,])*t(t(1-p)^(1-X[j,]))))*(-1/sqrt(sd2A)*(1-A^2/sd2A)))  #derivatives to sd(A) and sd(V)
    ret[j,2*nit+2]=sum(W*prod((t(t(p)^X[j,])*t(t(1-p)^(1-X[j,]))))*(-1/sqrt(sd2V)*(1-V^2/sd2V)))
  for(i in 1:nit){
    ret[j,i]= sum(W*derPai[,i]*X2[j,i]*prod((t(t(p)^X[j,])*t(t(1-p)^(1-X[j,])))[,-i]))
    ret[j,i+nit]= sum(W*derPvi[,i]*X2[j,i]*prod((t(t(p)^X[j,])*t(t(1-p)^(1-X[j,])))[,-i]))
  }
  }
}
return(ret)
}
