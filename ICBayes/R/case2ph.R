case2ph <-
function(L, 
R, 
status, 
xcov, 
x_user,
order,  
sig0, 
coef_range, 
a_eta,
b_eta,
knots,
grids,
niter){

#library(HI)

Ispline<-function(x,order,knots){
# M Spline function with order k=order+1. or I spline with order
# x is a row vector
# k is the order of I spline
# knots are a sequence of increasing points
# the number of free parameters in M spline is the length of knots plus 1.

k=order+1
m=length(knots)
n=m-2+k # number of parameters
t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m]) # newknots

yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
for (l in k:n){
    yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
}

yytem1=yy1
for (ii in 1:order){
   yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
   for (i in (k-ii):n){
      yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
   }
   yytem1=yytem2
}

index=rep(0,length(x))
for (i in 1:length(x)){
    index[i]=sum(t<=x[i])
}

yy=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))

if (order==1){
   for (i in 2:n){
      yy[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
   }
}else{
   for (j in 1:length(x)){
      for (i in 2:n){
         if (i<(index[j]-order+1)){
            yy[i-1,j]=1
         }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
            yy[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
         }else{
            yy[i-1,j]=0
         }
      }
   }
}
return(yy)
}

positivepoissonrnd<-function(lambda){ 
   samp = rpois(1, lambda)
   while (samp==0) {
      samp = rpois(1, lambda)
   }
   return(samp)
}

beta_fun <- function(x,j,beta,xx,te1,te2,sig0) 
{
  beta[j]<-x
  tt<-sum(xx[,j]*x*te1)-sum(exp(xx%*%beta)*te2)-x^2/sig0^2/2 # prior(beta)=N(0,sig0^2) 
  return(tt)     
}
ind_fun <- function(x,j,beta,xx,te1,te2,sig0) (x>-coef_range)*(x<coef_range)

### main routine ###
   L=matrix(L,ncol=1)
   R=matrix(R,ncol=1)
   R2=ifelse(is.na(R),0,R)
   xcov=as.matrix(xcov)
   p=ncol(xcov)
   status=matrix(status,ncol=1)
   n=nrow(L)
   u<-L*as.numeric(status==1)+R2*as.numeric(status==0)             
   v<-L*as.numeric(status==2)+R2*as.numeric(status==1)   
   obs=cbind(u,v)

   ## generate basis functions
   if (is.null(knots)) {knots<-seq(min(obs),max(obs),length=10)}
   k=length(knots)-2+order
   if (is.null(grids)) {grids<-seq(min(obs),max(obs),length=100)}
   kgrids=length(grids)
   G<-length(x_user)/p     # number of survival curves

   bisu=Ispline(t(obs[,1]),order,knots)
   bisv=Ispline(t(obs[,2]),order,knots)
   bgs=Ispline(grids,order,knots)

   eta=rgamma(1,a_eta,rate=b_eta)

   ## initial value
   gamcoef=matrix(rgamma(k, 1, 1),ncol=k)
   beta=matrix(rep(0,p),ncol=1)
   Lambdatu=t(gamcoef%*%bisu) # n x 1
   Lambdatv=t(gamcoef%*%bisv) # n x 1
   Lambdatg=t(gamcoef%*%bgs)

   parbeta=array(rep(0,niter*p),dim=c(niter,p))
   parsurv0=array(rep(0,niter*kgrids),dim=c(niter,kgrids))
   parsurv=array(rep(0,niter*kgrids*G),dim=c(niter,kgrids*G))
   pargam=array(rep(0,niter*k),dim=c(niter,k))
   pareta=array(rep(0,niter),dim=c(niter,1)) 
   parfinv=array(rep(0,niter*n),dim=c(niter,n))
   
   ## iteration
   iter=1
   while (iter<niter+1)
   { 
      # sample z, zz, w and ww
      z=array(rep(0,n),dim=c(n,1)); w=z
      zz=array(rep(0,n*k),dim=c(n,k)); ww=zz

      for (i in 1:n){
         if (status[i]==0){
            templam1=Lambdatu[i]*exp(xcov[i,]%*%beta)
            z[i]=positivepoissonrnd(templam1)
            zz[i,]=rmultinom(1,z[i],gamcoef*t(bisu[,i]))
         }else if (status[i]==1){
            templam1=(Lambdatv[i]-Lambdatu[i])*exp(xcov[i,]%*%beta)
            w[i]=positivepoissonrnd(templam1)
            ww[i,]=rmultinom(1,w[i],gamcoef*t(bisv[,i]-bisu[,i]))
         }
      }

      # sample beta
      te1=z*as.numeric(status==0)+w*as.numeric(status==1)
      te2=(Lambdatu*as.numeric(status==0)+Lambdatv*as.numeric(status==1)+Lambdatv*as.numeric(status==2))
      for (j in 1:p){
         beta[j]<-arms(beta[j],beta_fun,ind_fun,1,j=j,beta=beta,xx=xcov,te1=te1,te2=te2,sig0=sig0)
      }

      # sample gamcoef
      for (l in 1:k){
         tempa=1+sum(zz[,l]*as.numeric(status==0)*(bisu[l,]>0)+ww[,l]*as.numeric(status==1)*((bisv[l,]-bisu[l,])>0))
         tempb=eta+sum((bisu[l,]*as.numeric(status==0)+bisv[l,]*as.numeric(status==1)+bisv[l,]*as.numeric(status==2))*exp(xcov%*%beta)) 
         gamcoef[l]=rgamma(1,tempa,rate=tempb)
      }
      
      Lambdatu=t(gamcoef%*%bisu) # n x 1
      Lambdatv=t(gamcoef%*%bisv) # n x 1

      #sample eta
      eta=rgamma(1,a_eta+k, rate=b_eta+sum(gamcoef))

      pargam[iter,]=gamcoef
      parbeta[iter,]=beta
      pareta[iter]=eta
      ttt=gamcoef%*%bgs
      parsurv0[iter,]=exp(-ttt)
      if (is.null(x_user)){parsurv[iter,]=parsurv0[iter,]} else {
      A<-matrix(x_user,byrow=TRUE,ncol=p)
B<-exp(A%*%beta)
for (g in 1:G){
      parsurv[iter,((g-1)*kgrids+1):(g*kgrids)]=exp(-ttt*B[g,1])}
}

#calculate finv
Fu<-1-exp(-Lambdatu*exp(xcov%*%beta))   # n*1
Fv<-1-exp(-Lambdatv*exp(xcov%*%beta))   # n*1 
f_iter<-(Fu^(status==0))*((Fv-Fu)^(status==1))*((1-Fv)^(status==2)) # n*1, individual likelihood for each iteration
finv_iter<-1/f_iter            # n*1, inverse of individual likelihood for each iteration

      parfinv[iter,]=finv_iter

      iter=iter+1
   } # end iteration

   est<-list(parbeta=parbeta,
   parsurv0=parsurv0,
   parsurv=parsurv,
   parfinv=parfinv,
   grids=grids)
   est
}
