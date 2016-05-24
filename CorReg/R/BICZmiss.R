# ' BIC of Z with missing values in X
# ' 
# ' 


BICZmiss<-function(X=X,Z=Z,M=NULL,Bic_vide_vect=NULL,BicOld=NULL,methode=1,Zold=NULL,intercept=T,mixmod=NULL,nbit=100,nbclustmax=10,mixmodpenalty=T){
 
#   res=.Call( "BicZmissRcpp", X, Z, M, Bic_vide_vect, BicOld, methode, Zold, intercept, mixmod, nbit, nbclustmax, PACKAGE = "CorReg")
#   return(res)
# 
  if(is.null(M)){
    M=0*X
    M[is.na(M)]=1
    X[is.na(X)]=0
  }
  p=ncol(Z)
  if(is.null(Bic_vide_vect)){
    if(!is.null(mixmod)){
      Bic_vide_vect=rep(0,times=p)
      for(i in 1:p){
        Bic_vide_vect[i]=mixmod$BIC_vect[i]
      }
    }else{
       Bic_vide_vect=density_estimation(X=X)$BIC_vect
    }
  }
  if(is.null(Zold)| is.null(BicOld)){
    Zold=0*Z
    BicOld=Bic_vide_vect
  }
  BIC_vect=BicOld
  for(j in 1:p){#pour chaque variable (colonne)
    if(length(which(Z[,j]-Zold[,j]!=0))>0){#differences donc on recalcule Bicloc
      if(sum(Z[,j])==0){#colonne devenue vide
        BIC_vect[j]=Bic_vide_vect[j]
      }else{#Bic doit etre recalcule
        quiloc=(Z[,j]!=0)
        Xloc=as.matrix(X[,quiloc])
        Yloc=X[,j]
        Mloc=as.matrix(M[,quiloc])
        quibon=(M[,j]==0)
        Yloc=Yloc[quibon]
        Xloc=as.matrix(Xloc[quibon,])
        Mloc=as.matrix(Mloc[quibon,])
        mixmodloc=list()
        mixmodloc$nbclust=mixmod$nbclust[quiloc]
        mixmodloc$details=matrix(ncol=4,nrow=sum(mixmodloc$nbclust))
        k=1
        l=1
        m=quiloc[1]
        for(i in 1:nrow(mixmod$details)){
          if(Z[mixmod$details[i,4],j]!=0){#si la variable concernee est dans Xloc on recopie
            if(mixmod$details[i,4]==m){
              #on est sur la meme variable donc le numero ne change pas
            }else{#changement de variable donc il faut reindexer
              l=l+1
              m=mixmod$details[i,4]
            }
            mixmodloc$details[k,]=c(mixmod$details[i,1:3],l)
            k=k+1
          }
        }   
#         print(mixmodloc$details)
        res=OLS(X=Xloc,Y=Yloc,M=Mloc,intercept=intercept,methode=methode,miss=(sum(Mloc)>0),mixmod=mixmodloc,nbit=nbit,nbclustmax=nbclustmax,sigma=T)       
        BICLOC=c()  
        k=0#complexite du modele
        for (i in 1:nrow(Xloc)){#pour chaque point
          BICLOC=c(BICLOC,GM_Loglikelihood(Y=Yloc[i],X=Xloc[i,],M=Mloc[i,],B=res$beta,sigma=res$sigma,mixmod=mixmodloc,log=T,intercept=T))
          k=k+ncol(Xloc)+intercept+1#Beta+sigma dans tous les cas
          if(sum(Mloc[i,])>0){#trous donc on ajoute des parametres mixmod
            classesmank=Z[,j]*mixmodloc$nbclust*Mloc[i,] 
            #k=k+3sum(classmank)-nbfois ou classmank=1
            if(mixmodpenalty){
              k=k+prod(classesmank[classesmank>1]) #nombre de proportions estimees avec une par classe quand plusieurs classes          
              k=k+2*sum(classesmank) #moyennes et ecarts-types 
            }
          }            
        } 
        k=k/nrow(Xloc)
        BIC_vect[j]=-2*sum(BICLOC)+k*log(nrow(Xloc))
      }
    }
  }#fin de la boucle sur les colonnes    
  return(BIC_vect) 
}