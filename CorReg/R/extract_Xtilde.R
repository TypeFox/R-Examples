# ' extrait Xtilde
extract_Xtilde<-function(X=X,B=B){
  Z=B[-1,-1]
  Z=as.matrix(Z)
  Z=Z-diag(diag(Z))
  n=nrow(X)
  Z[Z!=0]=1
  quiI2=which(colSums(Z)!=0)
  quiI1=(1:ncol(Z))[-quiI2]
  X=as.matrix(X)
  return(X[,quiI2]-cbind(rep(1,times=n),X[,quiI1])%*%B[c(1,quiI1+1),quiI2])
}