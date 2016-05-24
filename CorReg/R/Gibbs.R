Gibbs<-function(last=FALSE,M=M,nbit=1,warm=0,mixmod=mixmod,X=X,comp_vect=comp_vect,missrow=missrow,quimiss=quimiss,
                Z=Z,Zc=Zc,alpha=alpha,sigma_IR=sigma_IR,nbclust_vect=nbclust_vect,Ir=Ir,loglik_bool=loglik_bool,Xout=FALSE,GibbsIR=TRUE){
   for(iter in 1:nbit){
      missrow_loc=1
      n=nrow(M)
      p=ncol(M)
      resmui=muiZ(p=ncol(Z),mixmod = mixmod,components=comp_vect[missrow_loc,],Z=Z,Zc=Zc,alpha=alpha,Ir=Ir,sigma_IR=sigma_IR)#on commence ligne 1
      mui=resmui$mui
      sigmai=resmui$sigmai
      Sigma=as.matrix(resmui$Sigma)
      loglik=rep(0,times=n)
      loglikfin=0
      nbmiss=sum(M)
      if(Xout){Xfin=X}
      for(i in 1:nbmiss){
         #          print(paste("i",i))
         miss=quimiss[i,]
         if(miss[1]!=missrow_loc){
            if(loglik_bool){#on a fini la ligne donc on calcule sa vraisemblance 
               loglik[missrow_loc]=loglikcond(X=X,mui=mui,Sigma=as.matrix(Sigma),M=M,i=missrow_loc,Zc=Zc)
            }
            missrow_loc=miss[1]#une fois la vraisemblance calculee on change de ligne
            resmui=muiZ(p=ncol(Z),mixmod = mixmod,components=comp_vect[missrow_loc,],Z=Z,Zc=Zc,alpha=alpha,Ir=Ir,sigma_IR=sigma_IR)#mise a jour
            mui=resmui$mui
            sigmai=resmui$sigmai
            Sigma=as.matrix(resmui$Sigma)
         }
         
         if(miss[2]%in%Ir){#redundant covariate missing
               X[miss[1],miss[2]]=Gibbs_X_ij_IR(X=X[miss[1],],sigma=sigma_IR[miss[2]],alpha=alpha[,miss[2]],GibbsIR=GibbsIR)
         }else{#missing right
            X[miss[1],miss[2]]=Gibbs_X_ij_IF(Z=Z,X=X[miss[1],],mui=mui,sigmai=sigmai,Sigma=Sigma,alpha=alpha,mixmod=mixmod,i=miss[1],j=miss[2],components=comp_vect[miss[1],])
         }
      }
      if(loglik_bool){#on a fini la derniere ligne donc on calcule sa vraisemblance 
         loglik[missrow_loc]=loglikcond(X=X,mui=mui,Sigma=as.matrix(Sigma),M=M,i=missrow_loc,Zc=Zc)
         for (i in 1:n){
            if(sum(M[i,])==0){#on veut la vraisemblance des lignes restantes (pleines)
               resmui=muiZ(p=ncol(Z),mixmod = mixmod,components=comp_vect[i,],Z=Z,Zc=Zc,alpha=alpha,Ir=Ir,sigma_IR=sigma_IR)#mise a jour
               mui=resmui$mui
               Sigma=as.matrix(resmui$Sigma)
               loglik[i]=loglikcond(X=X,mui=mui,Sigma=as.matrix(Sigma),M=M,i=i,Zc=Zc)
            }
         }
      }
      
      if(!last | (iter!=nbit)){
         #Imputation des classes####
         for(j in (1:p)[-Ir]){#pas besoin de classe a gauche car dans Gibbs on connait toujours tout le monde a droite (quitte a le supposer)
            #             print(paste("j",j))
            if(nbclust_vect[j]>1){
               for (i in 1:n){
                  #                   print("a")
                  vect=rmultinom(1, 1, tik(x=X[i,j],nbclust=nbclust_vect[j],mixmod=mixmod[[j]]))
                  #                   print("b")
                  #                   print(paste("vect",vect))
                  comp_vect[i,j]=match(1,vect)
               }
            }
         }
      } 
      loglikfin=sum(loglik)/nbit+loglikfin
      if(Xout){
         if (iter>1 & iter>(warm+1)){
            Xfin=Xfin+X/nbit
         }else if(iter==(warm+1)){
            Xfin=X/nbit
         }#sinon warm donc on ne fait rien
      }#optimiser en ne modifiant que les manquants  
   }
   #    print("c")
   #    print(X)
   #    print(comp_vect)
   if(loglik_bool){
      loglikfin=-2*loglikfin+missing_penalty(nbclust_vect=nbclust_vect,Z=Z,M=M,n=nrow(X),p=ncol(Z),Zc=Zc)#BIC adapte
   }
   result=list(X=X,comp_vect=comp_vect,loglik=loglikfin)    
   if(Xout){result$Xfin=Xfin
   }else{
      result$Xfin=X
   }
      return(result)
}