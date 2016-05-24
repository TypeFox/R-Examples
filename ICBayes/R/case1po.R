case1po <-
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

positivepoissonrnd<-function(lambda) 
{
    samp = rpois(1, lambda)
    while (samp==0) {
        samp = rpois(1, lambda)
    }
    return(samp)
}

beta_fun <- function(x,pb,xx,zz,lam, ksi,sig0) 
{
  p<-ncol(xx)
  tem<-matrix(c(pb,x),ncol=1)
  tt<-sum(zz*xx[,p]*x)-sum(lam*ksi*exp(xx%*%tem))-x^2/sig0^2/2 # prior(beta)=N(0,sig0^2) 
  return(tt)     
}
ind_fun <- function(x,pb,xx,zz,lam, ksi, sig0) (x>-coef_range)*(x<coef_range)



 ## obtain data
   T_obs=matrix(ifelse(status==0,L,R),ncol=1)
   status=matrix(status,ncol=1)
   n=nrow(T_obs)
   xcov=as.matrix(xcov)
   p=ncol(xcov)

   if (is.null(knots)) {knots<-seq(0,max(T_obs),length=10)}
   k=length(knots)-2+order
   if (is.null(grids)) {grids<-seq(0,max(T_obs),length=100)}
   kgrids=length(grids)
   G<-length(x_user)/p

   ## generate basis functions
   bis=Ispline(t(T_obs),order,knots) # used as basis functions in order to do parameter estimation
   bgs=Ispline(grids,order,knots)
   
   ## initial value
   eta=rgamma(1,a_eta,rate=b_eta)
   gamcoef=matrix(rgamma(k, 1, 1),ncol=k)
   beta=matrix(rep(0,p),ncol=1)
   Lambdat=t(gamcoef%*%bis) 
   ksi=array(-log(runif(n)), dim=c(n,1))

   parbeta=array(rep(0,niter*p),dim=c(niter,p))
   parsurv0=array(rep(0,niter*kgrids),dim=c(niter,kgrids))
   parsurv=array(rep(0,niter*kgrids*G),dim=c(niter,kgrids*G))
   pargam=array(rep(0,niter*k),dim=c(niter,k))
   pareta=array(rep(0,niter),dim=c(niter,1)) 
   parfinv=array(rep(0,niter*n),dim=c(niter,n))

   ## iterations
   iter=1
   while (iter<niter+1)
   {
      # sample z
      z=array(rep(0,n),dim=c(n,1))
      zz=array(rep(0,n*k),dim=c(n,k))
      for (i in 1:n){
         if (status[i]==1){
            templam=Lambdat[i]*exp(xcov[i,]%*%beta)*ksi[i]
            z[i]=positivepoissonrnd(templam)
            zz[i,]=rmultinom(1,z[i],gamcoef*t(bis[,i]))
         }
      }

      # sample ksi
      for (i in 1:n){
            ksi[i]=rgamma(1,1+z[i], rate=1+Lambdat[i]*exp(xcov[i,]%*%beta))
      }

      # sample gamcoef
      for (l in 1:k){
         tempa=1+sum(zz[,l]*(bis[l,]>0))
         tempb=eta+bis[l,]%*%(exp(xcov%*%beta)*ksi)
         gamcoef[l]=rgamma(1,tempa,rate=tempb)
      }
      Lambdat=t(gamcoef%*%bis) # n x 1

      #sample eta
      eta=rgamma(1,a_eta+k, rate=b_eta+sum(gamcoef))

      # sample beta
      for (l in 1:p){
         xx<-cbind(xcov[,-l],xcov[,l])
         pb<-beta[-l]
         beta[l]<-arms(beta[l],beta_fun,ind_fun,1,pb=pb,xx=xx,zz=z,lam=Lambdat,ksi=ksi, sig0=sig0)
      }

      pargam[iter,]=gamcoef
      parbeta[iter,]=beta
      pareta[iter]=eta
      parsurv0[iter,]=1/(1+gamcoef%*%bgs)
      if (is.null(x_user)){parsurv[iter,]=parsurv0[iter,]} else {
A<-matrix(x_user,byrow=TRUE,ncol=p)
B<-exp(A%*%beta)
for (g in 1:G){
parsurv[iter,((g-1)*kgrids+1):(g*kgrids)]=parsurv0[iter,]/(parsurv0[iter,]+(1-parsurv0[iter,])*B[g,1])}
}

#calculate finv
C<-Lambdat*exp(xcov%*%beta)
Ft<-C/(1+C)
f_iter<-(Ft^(status==1))*((1-Ft)^(status==0))  # n*1, individual likelihood for each iteration
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
