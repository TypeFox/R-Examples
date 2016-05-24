ubOSS <-
function(X, Y, verbose=TRUE){
  
  stopifnot(class(verbose) == "logical", all(unique(Y) %in% c(0, 1)))
  
  #only numeric features are allowed
  if(any(sapply(X,is.numeric)==FALSE))
    stop("only numeric features are allowed to compute nearest neighbors")
  
  S.X<-X
  S.Y<-Y
  i.1<-which(Y==1)
  N.1<-length(i.1)
  i.0<-which(Y==0)
  N.0<-length(i.0)
  if(N.1==0 | N.0==0) {
    cat("all instances of the same class \n")
    return(list(X=X,Y=Y))
  }
  
  #initially C contains all 1s from S and one random 0 obs
  id.C<-c(i.1,sample(i.0,1))	
  C.X<-X[id.C, ]
  C.Y<-Y[id.C]
  #use C to to build a 1-NN and classify all obs in S
  Y.knn<-knn(C.X, S.X, C.Y, k = 1)
  #move missclassified obs into C
  id.miss<-which(S.Y!=Y.knn)
  id.C<-c(id.C,id.miss)
  id.C <- sort(id.C)
  #id.C<-sample(id.C)
  C.X<-X[id.C, ]
  C.Y<-Y[id.C]
  #now C is consistent with S
  
  #remove from C 0s that are tomek links
  data<-ubTomek(C.X, C.Y, verbose)
  X<-data$X
  Y<-data$Y
  
  return(list(X=X,Y=Y))
}
