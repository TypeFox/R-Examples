case2probit <-
function(L,
R,
status,
xcov,
x_user,
order,
v0,
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


   L=matrix(L,ncol=1)
   R=matrix(R,ncol=1)
   R2=ifelse(is.na(R),0,R)
   status=matrix(status,ncol=1)
   n=nrow(L)
   xcov=as.matrix(xcov)
   p=ncol(xcov)
   err=1e-10
   u<-L*(status==1)+R2*(status==0)             
   v<-L*(status==2)+R2*(status==1)             
   obs=cbind(u,v)

   if (is.null(knots)) {knots<-seq(min(obs),max(obs),length=100)}
   if (is.null(grids)) {grids<-seq(min(obs),max(obs),length=10)}
   kgrids=length(grids)
   k=length(knots)-2+order
   G<-length(x_user)/p

   # initial values
   varbeta0=t(xcov)%*%xcov
   varbeta0inv=solve(varbeta0)
   Sigmabeta0=n*varbeta0inv
   invSigmabeta0=varbeta0/n
   intbeta0=matrix(rep(0,p),ncol=1)

   tt=obs[,1]*(status<=1)+obs[,2]*(status==2)
   bis=Ispline(t(tt),order,knots) # used as basis functions in order to do parameter estimation
   bisu=Ispline(t(obs[,1]),order,knots)
   bisv=Ispline(t(obs[,2]),order,knots)
   bisvy1=bisv[,status==1]
   bisuy1=bisu[,status==1] # for computational purpose
   bgs=Ispline(grids,order,knots)

   ## initial values
   eta=1
   m0=0
   gammapar=matrix(rep(1,k)/2,ncol=k); gamma0=-2;
   beta=matrix(rep(1,p),ncol=1) 
   Sigma0=diag(5,p) 
   invSigma0=solve(Sigma0)

   pareta=array(rep(0,niter),dim=c(niter,1))
   parsurv0=array(rep(0,niter*kgrids),dim=c(niter,kgrids))
   parsurv=array(rep(0,niter*kgrids*G),dim=c(niter,kgrids*G))
   pargamma0=array(rep(0,niter),dim=c(niter,1))
   pargamma=array(rep(0,niter*k),dim=c(niter,k))
   parbeta=array(rep(0,niter*p),dim=c(niter,p))
   parfinv=array(rep(0,niter*n),dim=c(niter,n))

   z=matrix(rep(0,n),ncol=1) # z=w_{i1} if y_i<=1 and w_{i2} if y_i=2. This definition matches the definition of t. 

   alphat=gamma0+t(gammapar%*%bis)
   alphau=gamma0+t(gammapar%*%bisu) # n x 1
   alphav=gamma0+t(gammapar%*%bisv) # n x 1
   alphag=gamma0+t(gammapar%*%bgs)

   ## iterations
   iter=1
   while  (iter<niter+1)
   {
      # sample z
      for (i in 1:n){
         if (status[i]==0){
            tempp1=pnorm(-alphau[i]-xcov[i,]%*%beta)
            temppp=min(1-err,runif(1)*(1-tempp1)+tempp1)
            z[i]=alphau[i]+xcov[i,]%*%beta+qnorm(temppp)
         }else if (status[i]==1){
            tempp1=pnorm(-alphau[i]-xcov[i,]%*%beta) 
            tempp2=pnorm(-alphav[i]-xcov[i,]%*%beta)
            temppp=min(1-err,runif(1)*(tempp1-tempp2)+tempp2)
            z[i]=alphau[i]+xcov[i,]%*%beta+qnorm(temppp)
         }else{
            tempp2=pnorm(-alphav[i]-xcov[i,]%*%beta)
            temppp=min(1-err,runif(1)*tempp2)
            z[i]=alphav[i]+xcov[i,]%*%beta+qnorm(temppp)
         }
      }

      # sample gamma0 
      tempw0=1/(v0+n)
      tempe0=tempw0*(v0*m0+ sum(z-t(gammapar%*%bis)-xcov%*%beta))
      gamma0=tempe0+rnorm(1)*sqrt(tempw0)

      zy1=z[status==1]; alphauy1=alphau[status==1]; alphavy1=alphav[status==1];

      # sample gamma's
      for (l in 1:k){
         tempindex=(bisvy1[l,]-bisuy1[l,]>0)
         if (sum(tempindex)>0){ #if sum(tempindex)=0, the likelihood does not contain gamma_l.
            tem1=gammapar[-l]%*%bisuy1[-l,tempindex>0]
            tem2=gammapar[-l]%*%bisvy1[-l,tempindex>0]
            tem3=matrix(bisvy1[l,tempindex>0]-bisuy1[l,tempindex>0],ncol=1)
            temperr=-(zy1[tempindex>0]-t(tem1)+t(tem2))/tem3
            temperq=max(temperr); cutoff=max(0, temperq);
            if (sum(bis[l,]^2)==0){
               gammapar[l]=cutoff -log(runif(1))/eta  # why add cutoff-? (paper p975)
            }else{
               tempb=1/(sum(bis[l,]^2))
               tempa=tempb*(bis[l,]%*%(z-gamma0-xcov%*%beta-t(gammapar[-l]%*%bis[-l,]))-eta)
               #tempf=pnorm(cutoff, tempa, sqrt(tempb))
               tempf=1-pnorm(cutoff, tempa, sqrt(tempb))
               tempu=runif(1); tempuu=min(1-tempf+tempu*tempf, 1-err);
               gammapar[l]=max(tempa+sqrt(tempb)*qnorm(tempuu),0)
            }
         }else{ 
            gammapar[l]=-log(runif(1))/eta
         }
      }

      alphat=gamma0+t(gammapar%*%bis)
      alphau=gamma0+t(gammapar%*%bisu) 
      alphav=gamma0+t(gammapar%*%bisv) 
      alphag=gamma0+t(gammapar%*%bgs)

      # sample beta
      tempsig=solve(invSigmabeta0+t(xcov)%*%xcov)
      tempbetam=tempsig%*%(invSigmabeta0%*%intbeta0+t(xcov)%*%(z-alphat))
      beta=t(chol(tempsig))%*%matrix(rnorm(p),ncol=1)+tempbetam

      # sample eta
      eta=rgamma(1, a_eta+k, rate=b_eta+sum(gammapar))

      # record parameters
      pareta[iter,]=eta
      pargamma0[iter,]=gamma0
      pargamma[iter,]=gammapar
      parbeta[iter,]=t(beta)
      parsurv0[iter,]=1-pnorm(alphag)
      if (is.null(x_user)){parsurv[iter,]=parsurv0[iter,]} else {
A<-matrix(x_user,byrow=TRUE,ncol=p)
B<-A%*%beta                 
for (g in 1:G){
      parsurv[iter,((g-1)*kgrids+1):(g*kgrids)]=1-pnorm(alphag+B[g,1])}
}

#calculate finv
Fu<-pnorm(alphau+xcov%*%beta)   # n*1
Fv<-pnorm(alphav+xcov%*%beta)   # n*1 
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
