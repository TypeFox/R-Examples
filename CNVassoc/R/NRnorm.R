NRnorm<-function(y,X,w,beta,sigma,variant,tol=10^-6,max.iter=1000,verbose=FALSE){

  if (missing(variant)) variant<-!apply(is.na(beta),1,any)
  beta<-t(apply(beta,1,function(x) if (any(is.na(x))) rep(x[1],length(x)) else x))

  J<-NCOL(w)
  n<-NROW(y)
  K<-NROW(beta)
  numbetes<-sum(ifelse(variant,J,1))

  if (max.iter<=0){
    derivatives<-hessianLinear(beta=beta,sigma=sigma,y=y,w=w,X=X,variant=variant)
    S<-derivatives$S
    H<-derivatives$H  
  }

  iter<-1
  error<-Inf

  while(error>tol & iter<max.iter){

    derivatives<-hessianLinear(beta=beta,sigma=sigma,y=y,w=w,X=X,variant=variant)
    S<-derivatives$S
    H<-derivatives$H

    beta.old<-beta
    sigma.old<-sigma
    beta.oldv<-matrix2vector(beta.old,variant)
    theta.oldv<-c(beta.oldv,sigma)
    thetav<-as.vector(theta.oldv-qr.solve(H)%*%S)
    betav<-thetav[1:numbetes]
    beta<-vector2matrix(betav,variant,J)
    sigma<-thetav[(numbetes+1):length(thetav)]

    error<-max(abs(thetav-theta.oldv))

    colnames(beta)<-paste("clust",1:J,sep="")
    rownames(beta)<-dimnames(X)[2][[1]]

    if (verbose){
      if (iter==1) cat("---- Newton-Raphson procedure ----\n")
      cat("Iter",iter,"Error",error,"\n")
      print(beta)
      cat("\n")
      cat("Sigma:\n")
      print(sigma)
    }

    iter<-iter+1

  }

  if (max.iter>0 & iter>max.iter) warning("Maximum iterations reached")

  return(list(beta=beta,sigma=sigma,variant=variant,score=S,hessian=H))

}
