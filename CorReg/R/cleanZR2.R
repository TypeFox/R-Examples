# ' To clean Z based on R2
# ' @export
# '@param Z binary adjacency matrix of the structure (size p)
# ' @param X the dataset
# ' @param R2min lower boundary for the structure (on R-squared value)
# ' @param methode parameter for OLS (matrix inversion) methode_BIC  parameter for OLS (matrix inversion) 1:householderQr, 2:colPivHouseholderQr
# ' @param adj boolean. Adjusted R-squared or classical one (if FALSE).
cleanZR2<-function(Z=Z,X=X,R2min=0.4,methode=1,adj=TRUE){
  p=ncol(Z)
  res=double(p)
  res_new=res
  quicol=which(colSums(Z)!=0)
  for(i in quicol){
    qui=which(Z[,i]!=0)
    ploc=length(qui)
    beta=OLS(X=X[,qui],Y=X[,i],intercept=T,methode=methode)$beta
    MSE=MSE_loc(Y=X[,i],X=as.matrix(X[,qui]),intercept=T,A=beta) #on met as.matrix pour les cas avec une seule colonne
    res[i]=1-(MSE)/(var((X[,i]))) 
    if(adj){
      res[i]=res[i]-(1-res[i])*ploc/(ncol(X)-ploc-1)
    }
    if(res[i]<R2min){#on supprime ce qui est pourri
      Z[,i]=0
    }else{
      res_new[i]=res[i]
    }
  }
  return(list(Z=Z,res_old=res,res_new=res_new))
}
