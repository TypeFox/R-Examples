Gibbs_miss<-function(X=X,alpha=alpha,nbit=1,Z=NULL,Zc=NULL,M=NULL,comp_vect=NULL,loglik_bool=FALSE,warm=0,Xout=FALSE,GibbsIR=TRUE,Ir=NULL,mixmod=NULL,sigma_IR=NULL){
   p=ncol(X)
   if(is.null(M)){
      M=0*X
      M[is.na(X)]=1
   }
   if(is.null(alpha)){
      alpha=hatB(Z = Z,X=X)
   }
   if(is.null(Z)){
      Z=matrix(0,ncol=p,nrow=p)
      Z[alpha[-1,]!=0]=1
   }
   if(is.null(Zc)){
      Zc=colSums(Z)
   }
   if(is.null(Ir)){
      Ir=which(Zc!=0)
   }
   if(is.null(sigma_IR)){
      sigma_IR=rep(0,times=p)#residus des regressions (stables pour gibbs)
      for (j in Ir){
         sigma_IR[j]=sd(X[,j]-X%*%alpha[-1,j])
      }
   }
   if(is.null(mixmod)){
      mixmod=density_estimation(X = X,detailed = TRUE)
   }
   missrow=sum(M)#nb d'individus a trou
   quimiss=which(M!=0,arr.ind = TRUE)
   if(is.null(comp_vect)){
      comp_vect=1+0*M
   }
   res=Gibbs(last=TRUE,M=M,nbit=nbit,warm=warm,mixmod=mixmod$details,X=X,comp_vect=comp_vect,missrow=missrow,quimiss=quimiss,
                       Z=Z,Zc=Zc,alpha=alpha,sigma_IR=sigma_IR,nbclust_vect=mixmod$nbclust,Ir=Ir,loglik_bool=loglik_bool,Xout=Xout,GibbsIR=GibbsIR)
   return(res)
}