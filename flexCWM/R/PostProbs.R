.PX_multi <- function(colXm,Xmod,z,nj,m,k,n){
  alpha <- NULL
  if(colXm>0){
    PXmult <- array(1,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
    for(u in 1:colXm)
      alpha[[u]] <- array((t(Xmod[[u]]) %*% z)/
                            matrix(rep(nj,m[u]),nrow=m[u],byrow=T),c(m[u],k),
                          dimnames=list(paste("level.",1:m[u],sep=""),paste("comp.",1:k,sep="")))
    names(alpha) <- names(Xmod)
    for(j in 1:k)
      for(i in 1:n){
        PXmult[i,j] <- 1
        for(u in 1:colXm)
          PXmult[i,j] <- PXmult[i,j]*dmultinom(Xmod[[u]][i,],size=1,prob=alpha[[u]][,j])
      }
  } else PXmult <- 1
  list(PX=PXmult, alpha=alpha)
}
.PX_norm <- function(colXn,Xnorm,modelXnorm,z,k,n,eps){
  if(colXn>0){
    PXnorm  <- array(1,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
    muX     <- array(0,c(colXn,k),dimnames=list(dimnames(Xnorm)[[2]],paste("comp.",1:k,sep="")))
    Sigma    <- array(0,c(colXn,colXn,k),dimnames=list(dimnames(Xnorm)[[2]],dimnames(Xnorm)[[2]],paste("comp.",1:k,sep="")))
    if(colXn>1){ #  More than one normal concomitant
      fitM <- mixture::m.step(data=Xnorm, covtype=modelXnorm, w=z, mtol=1e-10, mmax=10)
      for(j in 1:k){ 
        muX[,j]       <- fitM[[j]]$mu
        Sigma[,,j]    <- .fixSigma(fitM[[j]]$sigma,eps)
        log.det       <- determinant(as.matrix(Sigma[,,j]),logarithm=TRUE)$modulus
        PXnorm[,j]    <-exp(log((2*pi)^(-colXn/2)) + (-1/2)*log.det + (-1/2*mahalanobis(x=Xnorm, center=muX[,j], cov=Sigma[,,j], inverted=FALSE)))                                   
      }
    } else {
      fitM    <- mclust::mstep(modelName=substr(modelXnorm,1,1),data=Xnorm,z=z)
      muX[1,] <- fitM$parameters$mean
      if(substr(modelXnorm,1,1)=="E")
        Sigma[1,1,] <- rep(fitM$parameters$variance$sigmasq,k)
      if(substr(modelXnorm,1,1)=="V")
        Sigma[1,1,] <- fitM$parameters$variance$sigmasq
      for(j in 1:k)
        PXnorm[,j] <- dnorm(Xnorm,muX[1,j],sqrt(Sigma[1,1,j]))
    } 
  } else {
      PXnorm <- 1
      muX <- Sigma <- NULL
  }
  list(PX=PXnorm,mu=muX,Sigma=Sigma)
}
.fixSigma<- function(sigma,eps){
  es <- eigen(sigma)
  es$values[es$values<eps] <- eps
  es$vectors %*% diag(es$values) %*% t(es$vectors)
}
.PX_Poisson <- function(k,X,weights,n){
  if (!is.null(X)){
    lambdaX  <- array(0,c(ncol(X),k),dimnames=list(dimnames(X)[[2]],paste("comp.",1:k,sep="")))
    PX  <- array(1,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
    for (i in seq_len(ncol(X))){
      for(h in 1:k){
        lambdaX[i,h]   <- .PoissonX(X[,i],weights=weights[,h])$lambdaX
        PX[,h]  <- dpois(X[,i],lambda=lambdaX[i,h])              
      }
    }
      list(PX=PX, lambda=lambdaX)
  }else  list(PX=1, lambda=NULL)
  
}
.PX_bin <- function(k,X,weights,Xbtrials,n){
  if (!is.null(X)){
    p  <- array(0,c(ncol(X),k),dimnames=list(dimnames(X)[[2]],paste("comp.",1:k,sep="")))
    PX  <- array(1,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
    for (i in seq_len(ncol(X))){
      for(h in 1:k){
        p[i,h]    <- .BinomialX(X[,i],weights=weights[,h],m=Xbtrials[i])$p #stima p
        PX[,h]    <-  dbinom(X[,i],size=Xbtrials[i],prob= p[i,h])     
      }
    } 
    list(PX=PX, p=p)
  } else  list(PX=1, p=NULL)
 
}
.BinomialX <- function(X,weights,m){
  
  # X: a numerical vector for the unique covariate
  # weights: a numerical vector of weights the same length as X giving the weights to use for elements of X.
  
  # p.obs  <- weighted.mean(x=X/m,w=weights)
  
  f <- function(par,X,weights,m)
    -sum(weights*dbinom(x=X, size = m, prob=par, log = TRUE))
  
  res <- optimize(f=f, interval=c(0,1), X=X, weights=weights, m=m)
   
  return(list(loglik = -res$objective, p = res$minimum))  
  
}
.PoissonX <- function(X,weights){ 
  # X: a numerical vector for the unique covariate
  # weights: a numerical vector of weights the same length as X giving the weights to use for elements of X.
  lambda.obs  <- weighted.mean(x=X,w=weights)
  
  f <- function(par,X,weights){
    lambda <- par
    l <- -sum(weights*dpois(x=X, lambda = lambda, log = TRUE))
  }
  
  res <- optimize(f=f, interval=c(0,2*lambda.obs), X=X, weights=weights)
  
  loglik     <- -res$objective
  lambda.hat<-  res$minimum
  
  return(
    list(
      loglik  = loglik,
      lambdaX = lambda.hat
    )
  )  
  
}
