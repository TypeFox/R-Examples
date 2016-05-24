logLike.logistic<-function(param,y,X,w,variant){

  J<-NCOL(w)
 
  beta<-vector2matrix(betav=param,variant=variant,J=J)

  eta<-sapply(1:J,function(j){
      Xj<-cbind(X[,,j])
      betaj<-beta[,j,drop=FALSE]
      Xj%*%betaj
    })
    
  probs<-1/(1+exp(-eta))  
  ff<-sapply(1:J,function(j) dbinom(y,1,probs[,j]))
  gg<-rowSums(w*ff)
  return(sum(log(gg)))
}
