multComb<-function(n.vec){
  k<-length(n.vec)
  n<-sum(n.vec)

  our.combs<-combn(n,n.vec[1])

  next.step<-function(prev,n1,n){
    if(n1>1){
      new.cols<-t(combn((1:n)[-prev],n1))
    }
    if(n1==1){
      new.cols<-matrix((1:n)[-prev],ncol=1)
    }  
    t(cbind(matrix(prev,ncol=length(prev),nrow=max(1,nrow(new.cols)),byrow=T),new.cols))
  }

  for(i in 2:k){
    our.combs<-t(matrix(as.vector(apply(t(our.combs),1,next.step,n1=n.vec[i],n=n)),ncol=cumsum(n.vec)[i],byrow=T))
  }
  t(our.combs)
}