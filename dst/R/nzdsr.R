nzdsr<-function(x) 
{
  # Normalization of Dempster's rule of combination, Shafer method.
  # 
  # A normalization of the results from Dempster's rile of combination is done when there is a non zero mass allotted to the empty intersection. The method proposed by Shafer consist in dividing the results by 1 minus the mass allotted to the empty set.
  # 
  w12<-x$combination
  w1<-w12[,-1]
  mac<-w12[,1]
  i12<-x$I12
  nc=ncol(w12)-1  # tenir compte du cas où w12 a 1 ligne seul
  tri<-x$TRI
  con<-x$con
  ## remove empty set
  ind<-w12[tri[1],]
 # if ((sum(w12[1,1])!=0) & (sum(w12[1,-1])==0)) {
  if ((ind[1]!=0) & (sum(ind[-1])==0)) {
    ivide<-tri[1]  
    W2<-matrix(w1[tri[-1],],ncol=nc)
    ## calculate normalized masses
    MACC<-mac[tri[-1]]/(1-mac[ivide]) 
    MVIDE<-mac[ivide] 
  } 
  else {
    ivide<-0 
    W2<-matrix(w1,ncol=nc)
    MACC<-mac
    MVIDE<-0
  }
  W2<-cbind(matrix(MACC,ncol=1),W2)
  colnames(W2)<-colnames(w12)
  # W2 is the result of the normalization
  # the measure of conflict "con" is copied in the list
  list(DempsterRule=W2,con=con)
}