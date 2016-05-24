# ' calcul d'une structure initiale bas?e sur les corr?lations
# ' correl=cor(X_appr)
initial_cor3<-function(correl=correl,p2max=Inf,rmax=Inf){
  p=ncol(correl)
  Z=Matrix(0,nrow=p,ncol=p)
  qui1=which(matrix(runif(p^2),ncol=p)<correl)
  qui1=arrayInd(sample(qui1),c(p,p))#ordre al?atoire
  for(i in 1:nrow(qui1)){
    if(qui1[i,1]!=qui1[i,2] & (sum(Z[,qui1[i,1]])==0 & sum(Z[qui1[i,2],])==0) & sum(Z[,qui1[i,2]])<rmax & (length(which(colSums(Z)>0))<p2max | (length(which(colSums(Z)>0))==p2max & sum(Z[,qui1[i,2]])>0) ) ){
      Z[qui1[i,1],qui1[i,2]]=1
    }
  }
  return(Z)
}