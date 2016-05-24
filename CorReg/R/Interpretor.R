# 'another way to interprete the structure
Interpretor<-function (Z = Z, names = NULL, X = NULL,R2=T,R2_vect=NULL,adj=T) 
{
  if(R2){
    if(is.null(R2_vect)){
      R2_vect=R2Z(Z=Z,X=X,adj=adj,crit="R2")
    }
    ordre=order(R2_vect)
    X=X[,ordre]
    Z=Z[ordre,ordre]
    names=names[ordre]
  }
  
  if (is.null(names)) {
    if (!is.null(X)) {
      names = names(X)
    }else {
      print("no names")
      break
    }
  }
  
  groups = list()
  quiI2 = which(colSums(Z) != 0)
  compt = 1
  for (j in quiI2) {
    groups[[compt]] = c(names[j], names[which(Z[, j] != 0)])
    compt = compt + 1
  }
  if(R2){
    R2_vect=R2_vect[ordre]
    names(groups)=paste(1:length(groups),R2_vect[R2_vect!=0])
  }

  return(groups)
}