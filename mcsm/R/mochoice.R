mochoice=function(Niter=10^4){
#Model choice via MH steps
library(MASS)

#fertility data
y=log(as.vector(swiss[,1]))
X=as.matrix(swiss[,2:6])

#fast inverse matrix
inv=function(X){

  EV=eigen(X)
  EV$vector%*%diag(1/EV$values)%*%t(EV$vector)
  }

#subset of gam equal to 1
t1=function(gam){

   p=length(gam)
   if (sum(gam)!=0) order(gam)[(p-sum(gam)+1):p] else 0
   }
 
#marginal of model gam modulo constant
lpostw=function(gam,y,X,betatilde){

  n=length(y)
  qgam=sum(gam)
  t1gam=t1(gam)
  Xt1=cbind(rep(1,n),X[,t1gam])
  if (qgam!=0) P1=Xt1%*%inv(t(Xt1)%*%Xt1)%*%t(Xt1) else P1=matrix(0,n,n)
    -(qgam+1)/2*log(n+1)-n/2*log(t(y)%*%y-n/(n+1)*t(y)%*%P1%*%y-1/(n+1)*
    t(betatilde)%*%t(cbind(rep(1,n),X))%*%P1%*%cbind(rep(1,n),X)%*%betatilde)
  }

integerate=function(gamma){

  lga=length(gamma)
  vale=gamma[1]
  for (i in 2:lga)
    vale=vale+gamma[i]*2^(i-1)

  vale
  }

back2bin=function(k){

  gam=k%/%2^4;res=k%%2^4
  for (i in 3:0){
     logam=res%/%2^i;res=res%%2^i
     gam=c(logam,gam)
     }

  gam
  }

gocho=function(niter,y,X){

  lga=dim(X)[2];betatilde=lm(y~X)$coeff

  gamma=matrix(0,niter,lga)
  gamma[1,]=sample(c(0,1),lga,rep=T)
  for (i in 1:(niter-1)){

   j=sample(1:lga,1)
   gam0=gam1=gamma[i,];gam1[j]=1-gam0[j] 
   
   pr=lpostw(gam0,y,X,betatilde)
   pr=c(pr,lpostw(gam1,y,X,betatilde))
   pr=exp(pr-max(pr))

   gamma[i+1,]=gam0
   if (sample(c(0,1),1,prob=pr))
      gamma[i+1,]=gam1
   }

  gamma
  }

out=gocho(Niter,y,X)

# Best models
inde=apply(out,1,integerate)
cole=rep(0,length(unique(inde)))
for (i in 1:length(unique(inde))) cole[i]=(sum(inde==unique(inde)[i]))
besm=matrix(0,ncol=5,nrow=5)

for (i in 1:5)
  besm[i,]=back2bin(unique(inde)[order(-jitter(cole))[i]])

list(model=out,top=besm)
}
