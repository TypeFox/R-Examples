logLike.norm<-function(param,y,X,w,variant){

  J<-NCOL(w)
  K<-sum(ifelse(variant,J,1))

  beta<-vector2matrix(betav=param[1:K],variant=variant,J=J)
  sigma<-param[K+1]

  eta<-sapply(1:J,function(j){
      Xj<-cbind(X[,,j])
      betaj<-beta[,j,drop=FALSE]
      Xj%*%betaj
    })

  mu<-eta
  ff<-sapply(1:J,function(j) dnorm(y,mu[,j],sigma))
  gg<-rowSums(w*ff)
  return(sum(log(gg)))
  
}