dempster<-function(x,y) 
{
  ## x : tableau M lignes par K hypothèses + 1 colonne de masses
  ## y : tableau N lignes par K hypothèses + 1 colonne de masses
  ## -> inters : intersection des propositions
  ## -> combmasses : produit des vecteurs de masses
  ## -> doubles : remove duplicates of result of intersection
  nc<-ncol(x)-1
  N12<-inters(x,y)
  V12<-combmasses(x,y)
  ## transformer tableau : MxN lignes by K hypothesis
  N12<-aperm(N12,c(2,1,3))
  N12<-array(c(N12),c(dim(N12)[1],prod(dim(N12)[-1]))) ##tab k x (MxN)
  N12<-aperm(N12,c(2,1)) 
  ## remove duplicates of table
  W1<- doubles(N12)
  ## idendify les contributions to each subset
  I12<-dotprod(W1,aperm(N12,c(2,1)),g="&",f="==") 
  ## calculate total mass of each subset 
  MAC<-apply(I12*t(array(V12,dim(t(I12)))),1,sum)
  ## range hypothesis to check if the empty set is there
  TRI<-order(apply(W1,1,sum))
  ## Normalize masses if the empty set has a non zero mass
  z<- sum(W1[TRI[1],])
  ## identify the indice of the empty set
  if (z==0) IVIDE<-TRI[1] else IVIDE<-0 
  ## remove empty set
  if (z==0) W2<-matrix(W1[TRI[-1],],ncol=nc)  else W2<-matrix(W1,ncol=nc)
  ## remove the indices of the empty set
  if (z==0) K12<-I12[,TRI[-1]] else K12<-I12
  ## calculate normalized masses
  if (z==0) MACC<-MAC[TRI[-1]]/(1-MAC[IVIDE]) else MACC<-MAC
  ## calcul de MVIDE
  if (z==0) MVIDE<-MAC[IVIDE] else MVIDE<-0
  MVIDE<-array(MVIDE,c(nrow(W2),1))
  W2<-cbind(matrix(c(MVIDE,MACC),ncol=2),W2)
  return(W2)
}