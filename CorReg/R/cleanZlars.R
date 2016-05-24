# 'cleanZ by lars
# '
cleanZlars<-function(Z = Z, X = X, methode = 1, mode=c("MSE","BIC"),intercept=TRUE,K=NULL,groupe=NULL,select=c("lar", "lasso", "forward.stagewise", "stepwise",NULL)) 
{
  select=select[1]
  mode=mode[1]
  quicol = which(colSums(Z) != 0)
  for (i in quicol) {
    qui = which(Z[, i] != 0)
    Xloc = X[, qui]
    Yloc = X[, i]
    larsloc=lars(x=as.matrix(Xloc),y=Yloc,type=select,intercept=intercept)
    B_loc=meilleur_lars(lars=larsloc,X=Xloc,Y=Yloc,mode=mode,intercept=intercept,K=K,groupe=groupe)$A
    Z[qui,i][B_loc[-1]==0]=0#on repercute les 0
  }
  return(Z)
}