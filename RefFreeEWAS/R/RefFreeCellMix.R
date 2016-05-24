projectMix <- function(Y, Xmat, nonnegative=TRUE, sumLessThanOne=TRUE, lessThanOne=!sumLessThanOne){

  nCol = dim(Xmat)[2]
  nSubj = dim(Y)[2]

  mixCoef = matrix(0, nSubj, nCol)
  rownames(mixCoef) = colnames(Y)
  colnames(mixCoef) = colnames(Xmat)

  if(nonnegative){

    if(sumLessThanOne){
      Amat = cbind(rep(-1,nCol), diag(nCol))
      b0vec = c(-1,rep(0,nCol))
    }
    else if(lessThanOne){
      Amat = cbind(-diag(nCol), diag(nCol))
      b0vec = c(rep(-1,nCol),rep(0,nCol))
    }
    else{
      Amat = diag(nCol)
      b0vec = rep(0,nCol)
    }

    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i]))
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec)$sol
    }
  }
  else{
    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i]))
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
    }
  }

  return(mixCoef)
}


RefFreeCellMix <- function(Y,mu0=NULL,K=NULL,iters=10,Yfinal=NULL,verbose=TRUE){
  if(is.null(mu0)){
    if(K==1) {
       if(!is.null(Yfinal)) Y <- Yfinal
       n <- dim(Y)[2]

       mu <- matrix(apply(Y,1,mean,na.rm=TRUE),ncol=1)
       omega <- matrix(1, n, 1)
       o <- list(Mu=mu, Omega=omega)
       class(o) <- "RefFreeCellMix"
       return(o)
    }
    else mu0 <- RefFreeCellMixInitialize(Y,K=K,method="ward")
  }
  for(i in 1:iters){
    flag <- !apply(is.na(mu0),1,any)
    omega <- projectMix(Y[flag,],mu0[flag,])
    mu <- projectMix(t(Y),omega,sumLessThanOne=FALSE)
    if(verbose) print(summary(abs(as.vector(mu-mu0))))
    mu0 <- mu
  }
  if(!is.null(Yfinal)){
    mu <- projectMix(t(Yfinal),omega,sumLessThanOne=FALSE)
  }

  o <- list(Mu=mu, Omega=omega)
  class(o) <- "RefFreeCellMix"
  o
}

print.RefFreeCellMix <- function(x,...){
  cat("Reference Free Deconvolution\n\n")
  cat("Mu: ",dim(x$Mu)[1]," cpgs x ",dim(x$Mu)[2],"cell types\n")
  cat("Omega :",dim(x$Omega)[1]," subjects x ",dim(x$Omega)[2],"cell types\n")
}

summary.RefFreeCellMix <- function(object,...){
  list(Mu=apply(object$Mu,2,summary), Omega=apply(object$Omega,2,summary), MuCorr=cor(object$Mu))
}

RefFreeCellMixArray <- function(Y,Klist=1:5,iters=10,Yfinal=NULL,verbose=FALSE,dist.method = "euclidean",...){
   D <- dist(t(Y),method=dist.method)
   hc <- hclust(D,...)

   rfcmArray <- list()
   nK <- length(Klist)
   for(r in 1:nK){
     cat("Fitting K =",Klist[r],"\n")
     if(Klist[r]==1){
       rfcmArray[[r]] <- RefFreeCellMix(Y,K=1,Yfinal=Yfinal)
     }
     else{
       rfcmArray[[r]] <- RefFreeCellMix(Y,mu0=RefFreeCellMixInitialize(Y,K=Klist[r],Y.Cluster=hc),
         Yfinal=Yfinal, verbose=verbose)
     }
   }
   names(rfcmArray) <- Klist
   rfcmArray
}

deviance.RefFreeCellMix <- function(object, Y, Y.oob=NULL, EPSILON=1E-9, 
  bootstrapIterations=0, bootstrapIndices=NULL, ...){

  N <- dim(Y)[2]
  if(bootstrapIterations>0){# Do the bootstrap and replace x (but initialize with x$Mu)
    if(is.null(bootstrapIndices)){
      boots <- sample(1:N, N, replace=TRUE)
    }  
    else {
      boots <- bootstrapIndices
    }
    Y.oob <- Y[,-unique(boots)]
    Y <- Y[,boots]
    if(dim(object$Mu)[2]==1) {
      object <- RefFreeCellMix(Y,K=1,iters=bootstrapIterations,verbose=FALSE)
    }
    else{
      object <- RefFreeCellMix(Y,mu0=object$Mu,iters=bootstrapIterations,verbose=FALSE)
    }
  }

  Y.mu <- object$Mu %*% t(object$Omega)
  R <- Y-Y.mu
  Y.n <- apply(!is.na(Y),1,sum)
  Y.SSQ <- apply(R*R,1,sum,na.rm=TRUE)
  logSigma2 <- log(pmax(EPSILON,Y.SSQ)) - log(Y.n)

  if(!is.null(Y.oob)){
    Omega.oob <- projectMix(Y.oob,object$Mu)
    Y.mu <- object$Mu %*% t(Omega.oob)
    R.oob <- Y.oob-Y.mu
    n.oob <- apply(!is.na(Y.oob),1,sum)
    SSQ.oob <- apply(R.oob*R.oob,1,sum,na.rm=TRUE)
    N <- dim(Y.oob)[2]
  }
  else{
    SSQ.oob <- Y.SSQ
    n.oob <- Y.n
  }

  sum( n.oob*log(2*pi)+n.oob*logSigma2+SSQ.oob/exp(logSigma2))/N
}

RefFreeCellMixInitialize <- function(Y,K=2,Y.Distance=NULL, Y.Cluster=NULL, 
    largeOK=FALSE, dist.method = "euclidean", ...){

   if(!is.matrix(Y) | !is.numeric(Y)){
        stop("Y is not a numeric matrix\n")
   }
   n <- dim(Y)[2]

   if(is.null(Y.Cluster)){
     if(is.null(Y.Distance)){
        if(n>2500 & !largeOK){
           stop("Y has a large number of subjects!  If this is what you really want, change 'largeOK' to TRUE\n")
        }
        Y.Distance <- dist(t(Y),method=dist.method)
      }
      Y.Cluster <- hclust(Y.Distance,...)
   } 

   classes <- cutree(Y.Cluster, K)
   s <- split(1:n,classes)

   sapply(s, function(u) apply(Y[,u,drop=FALSE],1,mean,na.rm=TRUE))
}

RefFreeCellMixArrayDevianceBoot <- function(rfArray, Y, EPSILON=1E-9, bootstrapIterations=5){
  N <- dim(Y)[2]
  boots <- sample(1:N, N, replace=TRUE)
  sapply(rfArray, deviance.RefFreeCellMix, Y=Y, EPSILON=EPSILON, 
    bootstrapIterations=bootstrapIterations,bootstrapIndices=boots)
} 

RefFreeCellMixArrayDevianceBoots <- function(rfArray, Y, R=5, EPSILON=1E-9, bootstrapIterations=5){
  dv <- sapply(rfArray, deviance, Y=Y)
  nK <- length(dv)
  devs <- matrix(NA,R,nK)
  for(r in 1:R){
     if(r %% 10==0) cat("Bootstrap",r,"\n")
     devs[r,] <- RefFreeCellMixArrayDevianceBoot(rfArray, Y, 
       EPSILON=EPSILON, bootstrapIterations=bootstrapIterations)
  }
  out <- rbind(dv,devs)
  rownames(out) <- 0:R
  out
} 

