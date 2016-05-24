EMlogistic<-function(y,X,w,beta,tol=10^-6,max.iter=1000,verbose=FALSE){

  variant<-!apply(is.na(beta),1,any)
  beta<-t(apply(beta,1,function(x) if (any(is.na(x))) rep(x[1],length(x)) else x))

  P<-ncol(X)
  J<-ncol(w)
  N<-length(y)

  error<-Inf
  iter<-1

  latent<-rep(1:J,each=N)

  X.expand<-NULL
  for (p in 1:P){
    xp<-cbind(as.vector(X[,p,]))
    colnames(xp)<-paste("V",p,sep="")
    if (variant[p]){
      xp<-model.matrix(~xp:as.factor(latent)-1)
      colnames(xp)<-paste(rep("V",J),p,sep="")
    }
    X.expand<-cbind(X.expand,xp)
  }

  y.expand<-rep(y,J)

  while(error>tol & iter<max.iter){

    beta.old<-beta

    eta<-sapply(1:J,function(j){
        Xj<-cbind(X[,,j])
        betaj<-beta[,j,drop=FALSE]
        Xj%*%betaj
      })

    p<-1/(1+exp(-eta))


    dY<-sapply(1:J, function(j) dbinom(y,1,p[,j]))
    wp<-dY*w
    wp<-wp/rowSums(wp)

    coeffic<-logistic.fit(x=X.expand,y=y.expand,weights=as.vector(wp))

    beta<-vector2matrix(betav=coeffic,variant=variant,J=J)

    error<-max(abs(beta-beta.old))

    beta.old<-beta

    colnames(beta)<-paste("clust",1:J,sep="")
    rownames(beta)<-dimnames(X)[2][[1]]

    if (verbose){
      if (iter==1) cat("--- EM procedure ---\n")
      cat("Iter",iter,"Error",error,"\n")
      print(beta)
      cat("\n")
    }

    iter<-iter+1

  }

  return(list(beta=beta,variant=variant))

}

