ubCNN <-
function(X,Y,k=1,verbose=T){
  
  #only numeric features are allowed
  is.not.num<-which(sapply(X,is.numeric)==FALSE)
  if(length(is.not.num)>0)
    stop("only numeric features are allowed to compute nearest neighbors")
  
  S.X<-X
  S.Y<-Y
  i.1<-which(Y==1)
  i.0<-which(Y==0)
  N.1<-length(i.1)
  N.0<-length(i.0)
  if(N.1==0 | N.0==0) {
    if(verbose) cat("All instances of the same class \n")
    return(list(X=X,Y=Y))
  }
  
  #initially C contains all 1s from S and one random 0 obs
  id.C<-c(i.1,sample(i.0,1))	
  C.X<-X[id.C,]
  C.Y<-Y[id.C]
  #use C to to build a 1-NN and classify all obs in S
  Y.knn<-knn(C.X,S.X,C.Y,k)
  #move missclassified obs into C
  id.miss<-which(S.Y!=Y.knn)
  id.C<-c(id.C,id.miss)
  # id.C<-sample(id.C)
  X<-X[id.C, ]
  Y<-Y[id.C]
  #now C is consistent with S
  
  return(list(X=X,Y=Y))
}
