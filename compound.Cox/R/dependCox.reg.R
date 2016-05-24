dependCox.reg <-
function(t.vec,d.vec,X.vec,alpha,var=TRUE){
  
alpha=max(alpha,0.01) ### prevent unstability ###

n=length(t.vec)
t.ot=sort(t.vec)
d1.ot=(d.vec)[order(t.vec)]
d2.ot=(1-d.vec)[order(t.vec)]
x.ot=X.vec[order(t.vec)]
n1=sum(d1.ot)
n2=sum(d2.ot)

l.func=function(dL){
  b1=dL[n+1];b2=dL[n+2]
  dR1=dR2=rep(0,n)
  dR1[d1.ot==1]=exp(dL[1:n1])
  dR2[d2.ot==1]=exp(dL[(n1+1):n])
  R1=cumsum(dR1)
  R2=cumsum(dR2)
  E1=exp(-R1*exp(b1*x.ot))^(-alpha)
  E2=exp(-R2*exp(b2*x.ot))^(-alpha)
  eta1=E1/(E1+E2-1)
  eta2=E2/(E1+E2-1)
  l=sum( (b1*x.ot+log(eta1)+log(dR1))[dR1>0] )+sum( (b2*x.ot+log(eta2)+log(dR2))[dR2>0] )
  l=l-sum( log(E1+E2-1)/alpha )
  -l  
}

dL0=c(log( c(rep(1,n1)/n,rep(1,n2)/n) ),0,0)

if(var==FALSE){
  res=nlm(l.func,p=dL0,hessian=var)
  c(beta_hat=res$estimate[n+1])
}
else{
  res=nlm(l.func,p=dL0,hessian=var)
  beta_d=res$estimate[n+1]
  SD=sqrt(solve(res$hessian)[n+1,n+1])
  Z_d=beta_d/SD
  P_d=1-pchisq(Z_d^2,df=1)
  c(beta_hat=beta_d,SD=SD,Z=Z_d,P=P_d)
}

}
