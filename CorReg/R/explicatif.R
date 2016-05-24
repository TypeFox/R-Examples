# ' calcul des modeles explicatifs
explicatif<-function(X=X,Y=Y,Z=Z,type=c("lasso", "lar", "forward.stagewise", "stepwise"),intercept=TRUE){
  X=X[,colSums(Z)==0]
  res=lars(x=X,y=Y,type=type,intercept=intercept)
  return(res)
}