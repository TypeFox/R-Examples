# 'On estime les coefficients du meilleur mod?le obtenu (meilleur au sens du crit?re choisi)
# ' 
# ' 
meilleur_lars<-function(lars=lars,X=X,Y=Y,mode=c("MSE","BIC"),intercept=TRUE,K=NULL,groupe=NULL,Amax=NULL,OLS=TRUE){
  mode=mode[1]
  X=as.matrix(X)
  if(is.null(Amax)){Amax=ncol(X)+1}
  if(mode=="BIC"){
    vrais_vect=c()
    coef=as.matrix(coef(lars))
    coefbin=coef
    coefbin[coef!=0]=1
    coef=as.matrix(coef[rowSums(coefbin)<=Amax,])
    A=rep(0,times=ncol(X)+intercept)
    if(!is.null(nrow(coef))){
      for(i in 1:nrow(coef)){
        qui=which(coef[i,]!=0)
        vrais_vect[i]=OLS(X=X[,qui],Y=Y,intercept=intercept,Bic=T)$BIC
      }
      qui=which.min(vrais_vect)
      BIC=vrais_vect[qui]
      qui=which(coef[qui,]!=0)
      A[c(intercept,qui+intercept)]=OLS(X=X[,qui],Y=Y,intercept=intercept)$beta
    }
    return(list(A=A,BIC=BIC))    
  }else{#on regarde les MSE et on prend le meilleur
    K=min(K,nrow(X))
    vrais_vect=c()
    coef=as.matrix(coef(lars))
    coefbin=coef
    coefbin[coef!=0]=1
    coef=as.matrix(coef[rowSums(coefbin)<=Amax,])
    A=rep(0,times=ncol(X)+intercept)
    MSE=var(Y)
    if(!is.null(nrow(coef))){
      for(i in 1:nrow(coef)){
        qui=which(coef[i,]!=0)
        vrais_vect[i]=CVMSE(X=X[,qui],intercept=intercept,Y=Y,K=K,groupe=groupe)
      }
      qui=which.min(vrais_vect)
      MSE=vrais_vect[qui]
      qui=which(coef[qui,]!=0)
      if(length(qui)>0){
         A[c(intercept,qui+intercept)]=OLS(X=as.matrix(X[,qui]),Y=Y,intercept=intercept)$beta
      }#sinon ne garde rien donc on laisse A sur sa valeur initiale nulle
   }
    return(list(A=A,CVMSE=MSE))    
  }
}