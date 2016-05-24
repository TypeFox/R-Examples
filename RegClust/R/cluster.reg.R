cluster.reg <-function(Y,X,loop=1000){
X<-data.frame(X)
Y<-data.frame(Y)

if(ncol(X) > 3)stop(paste("too many covariates"))
if(nrow(X) != nrow(Y))stop(paste("unequal number of observations"))
if(ncol(X) >= nrow(X))stop(paste("too many covariates: number of covariates larger than sample size"))

start<-0
Numcov<-ncol(X)
NumCol<-seq(1:ncol(Y))
x<-seq(1,nrow(X))
size<-dim(X)[1]
bootstrap<-1
NoClust<-ncol(X)
loop<-loop
probSamp<-rep(1/nrow(X),nrow(X))

cluster<-matrix(0,nrow=bootstrap,ncol=length(NumCol))
ColParm<-(Numcov+1)+1+1
parm<-matrix(0,nrow=bootstrap*NoClust,ncol=ColParm)
dim(parm)
samp<-1

  #define the design matrix indep#
   
   intercept<-rep(1,nrow(X))
   indep<-as.matrix(data.frame(intercept,X))
   is.numeric(indep)
   
  #Initializing #
   
   beta<-matrix(rnorm(NoClust*(Numcov+1),0,0.1),byrow=T,ncol=Numcov+1)
   s<-exp(rnorm(NoClust,0,0.1))
   p<-runif(NoClust)
   sumP<-sum(p)
   p<-p/sumP

   lhood=c()

 
  #EM algorithm for each subset#

   u<-matrix(NA,nrow=nrow(beta),ncol=length(NumCol))
   for (w in 1:loop){
     numerator=matrix(0,nrow(beta),length(NumCol))
     for (g in 1:(length(NumCol))){
       for (k in 1:nrow(beta)){
   numerator[k,g]=log(p[k])-nrow(Y)*log(s[k])-sum((Y[,g+start]-indep%*%beta[k,])^2/(2*s[k]^2))

       }#end for (k in ...)
 maxNum=max(numerator[,g])
       total=0
       for (k in 1:nrow(beta)){
         total=total+exp(numerator[k,g]-maxNum)
       }
       for (k in 1:nrow(beta)){
       u[k,g]=exp(numerator[k,g]-maxNum-log(total))
 }
     } #end for (g in ...)


  ###likelihood in E-step###

   aa<-matrix(0,nrow(u),ncol(u))
   for (g in 1:ncol(u)){
     for (k in 1:nrow(u)){
       aa[k,g]<-log(p[k])*u[k,g]
     }
   }
   l1=sum(aa)
 
   L2=c()
   L21=matrix(0,nrow(u),ncol(u))
   for (g in 1:ncol(u)){
     for (k in 1:nrow(u)){
  for (i in 1:nrow(Y)){
         #L2[i]=log(s[k]^2)+(eye[i,g+start]-indep[i,]%*%(beta[k,]))^2/s[k]^2
         L2[i]=log(s[k]^2)+(Y[i,g+start]-indep[i,]%*%(beta[k,]))^2/s[k]^2
       }
       L21[k,g]=sum(L2)
     }
   }
   l2=sum(u*L21)
   l<-l1-0.5*l2
   lhood<-c(lhood,l) 

 ###M-step###
 #update pi (p)#
   p<-matrix(0,1,nrow(u))
   for (k in 1:nrow(u)){
   if (sum(u[k,])==0)  u[k,]<-10^(-40)
   p[k]<-sum(u[k,])/ncol(u)
   }

 #update beta#
   beta=matrix(rep(0,nrow(beta)*ncol(beta)),nrow=nrow(beta))
   for (k in 1:nrow(beta)){
     beta1=matrix(rep(0,ncol(beta)*ncol(beta)),nrow=ncol(beta))
     beta2=rep(0,ncol(beta))
     for (g in 1:ncol(u)){
       beta1<-u[k,g]*t(indep)%*%indep+beta1
       beta2<-u[k,g]*t(indep)%*%Y[,g+start]+beta2
     }
     beta[k,]<-solve(beta1)%*%beta2
   }

 #update sigma#
   sigma1=c()
   tempsigma=c()
   for (k in 1:nrow(u)){
     for (g in 1:ncol(u)){
       tempsigma[g]<-u[k,g]*t(Y[,g+start]-indep%*%(beta[k,]))%*%(Y[,g+start]-indep%*%(beta[k,]))
     }
     sigma1[k]=sum(tempsigma)
   }

   sigma2=c()
   for (k in 1:nrow(beta)){
     sigma2[k]=nrow(Y)*sum(u[k,])

   }

   sigma=c()
   for (k in 1:nrow(beta)){
     sigma[k]=sigma1[k]/sigma2[k]
   }
   s<-sqrt(sigma)

   if (w==1) {
     l[w]=l
   }
   else if (w >=2){
     #cat("lh ",w," ",lhood[w]," w-1 ",lhood[w-1],"\n")
     if ((lhood[w]-lhood[w-1])<10^(-40)) {  
       cat("Converged at iteration ",w,", BIC=",2*lhood[w]-Numcov*log(nrow(Y)),"\n")
 break
}
   }

if(w == loop){  cat("Did not converge, BIC=",2*lhood[w]-Numcov*log(nrow(Y)),"\n")}       

   w<-w+1
   }#end for (w in 1:1000)

 #output parameters#
 if (samp==1) {
  parm<-cbind(beta,t(rbind(sigma,p)))
 }

 #identify clusters that each Zernike values is in#
   for (g in 1:ncol(u)){
     max<-u[1,g]
     index<-1
     for (k in 2:nrow(u)){
       if (u[k,g]>max){
         max<-u[k,g]
         index<-k
       }#end if (u[k,g]>max)
     }#end for (k in 2:nrow(u))
     cluster[samp,g]<-index
   }#end for (g in 1:ncol(u))
   subClust<-matrix(0,nrow=nrow(u),ncol=ncol(u))
   for (k in 1:nrow(u)){
     m<-0
     for (g in 1:ncol(u)){
       if (cluster[samp,g]==k) {       
         m<-m+1
         subClust[k,m]<-g
         }#end if (cluster[samp,g]...)
     }#end for (g in 1:ncol...)
     position<-c()
     if (m>=1) {
       for (kk in 1:m){
         position[kk]<-subClust[k,kk]+Numcov+1
       }#end for (kk in 1:m)
     }#end if (m>=1)
   }#end for (k in ...)

clusterIndex<-unique(cluster[1,])
for (j in 1:length(clusterIndex))
{jj<-0
r2<-0
for (i in 1:length(NumCol))
{
if (cluster[1,i]==clusterIndex[j])
{
jj<-jj+1
rr<-sum((indep%*%parm[clusterIndex[j],1:ncol(beta)]-mean(Y[,(i+start)]))^2)
r0<-sum((Y[,(i+start)]-mean(Y[,(i+start)]))^2)
r2<-r2+rr/r0
}
}
}

colnames(cluster)<-colnames(Y)
colnames(parm)<-c("Intercept",colnames(X),"sigma","pi")
rownames(parm)<-paste("cluster ",sort(unique(cluster[1,])),sep='')

output=list(cluster=cluster,param=parm,likelihood=lhood[w],BIC=2*lhood[w]-Numcov*log(nrow(Y)))

output
}