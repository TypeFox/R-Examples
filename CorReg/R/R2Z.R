# ' Estimates R2 of each subregression
# '@param Z binary adjacency matrix of the structure (size p)
# ' @param X the dataset
# ' @param methode parameter for OLS (matrix inversion) methode_BIC  parameter for OLS (matrix inversion) 1:householderQr, 2:colPivHouseholderQr
# ' @param adj boolean to choose between adjusted R-squared and classical one
# ' @param crit to choose between the R-squared and the F statistic (p-value)
# ' @export 
R2Z<-function(Z=Z,X=X,methode=1,adj=F,crit=c("R2","F","sigmaX")){
  p=ncol(Z)
  res=rep(0,times=p)
  crit=crit[1]
  quicol=which(colSums(Z)!=0)
  for(i in quicol){
    qui=which(Z[,i]!=0)
#     ploc=length(qui)
#     beta=OLS(X=X[,qui],Y=X[,i],intercept=T,methode=methode)$beta
#     MSE=MSE_loc(Y=X[,i],X=as.matrix(X[,qui]),intercept=T,A=beta) #on met as.matrix pour les cas avec une seule colonne
#     res[i]=1-(MSE)/(var((X[,i])))  
    Xloc = X[, qui]
    Yloc = X[, i]
    lmloc=lm(Yloc~.,data=data.frame(Xloc))
    summar=summary(lmloc)
    if(crit=="R2"){
      if(adj){
  #       res[i]=res[i]-(1-res[i])*ploc/(ncol(X)-ploc-1)  
        res[i]=as.numeric(summar[9])
      }else{
        res[i]=as.numeric(summar[8])
      }
    }else if(crit=="F"){#p-value du test F global de Fisher
      res[i]= pf(summar$fstatistic[1],summar$fstatistic[2],summar$fstatistic[3],lower.tail=FALSE)   
    }else{#sigmaX
       res[i]=sd(summar[[3]])
    }
  }
  return(res)
}
