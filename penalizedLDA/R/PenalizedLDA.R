PenalizedLDA <-
function(x, y, xte=NULL,  type="standard",  lambda,  K=2, chrom=NULL, lambda2=NULL, standardized=FALSE, wcsd.x=NULL, ymat=NULL, maxiter=20, trace=FALSE){
  if(sum(1:length(unique(y)) != sort(unique(y)))>0) stop("y must be a numeric vector, with values as follows: 1, 2, ....")
  if(sum(is.na(x))>0 || sum(is.na(y))>0 || (!is.null(xte) && sum(is.na(xte))>0)) stop("No missing values allowed!!!")
  if(K>=length(unique(y))) stop("Can have at most K-1 components of K unique classes")
  yclass <- y
  y <- ymat
  if(is.null(ymat)) y <- MakeYMat(yclass)
  if(type=="ordered" && is.null(lambda2)) stop("For type 'ordered', lambda2 must be specified.")
  xorig <- x
  
  if(!standardized){
    if(is.null(wcsd.x)){
      if(length(y)<=200){
        wcsd.x <- wcsd.matrix(x,y)#apply(x, 2, wcsd, yclass)
      } else {
        wcsd.x <- apply(x,2,wcsd,yclass)
      }
      if(min(wcsd.x)==0) stop("Some features have 0 within-class standard deviation.")
    }
    if(!is.null(xte)) xte <- scale(xte, center=apply(x,2,mean), scale=wcsd.x)
    x <- scale(x, T, scale=wcsd.x)
  }
  sqrt.sigma.bet <- t(scale(y, F, sqrt(apply(y, 2, sum))))%*%x/sqrt(nrow(x))
  while(sum(is.na(sqrt.sigma.bet))>0){
    sqrt.sigma.bet <- t(scale(y, F, sqrt(apply(y, 2, sum))))%*%x/sqrt(nrow(x))
    cat("retrying", fill=TRUE)
  }
  penpca <- PenalizedPCA(x=sqrt.sigma.bet, lambda=lambda,  K=K, type=type, chrom=chrom, lambda2=lambda2, maxiter=maxiter, trace=trace)
  Usparse <- penpca$v#matrix(penpca$v, ncol=K)
  if(K==1) Usparse <- matrix(Usparse,ncol=1)
  if(sum(is.na(Usparse))>0){
    Usparse[is.na(Usparse)] <- 0
    #Usparse <- matrix(Usparse,ncol=K)
    if(K==1) Usparse <- matrix(Usparse,ncol=1)
  }
  xtranssparse <- x%*%Usparse# matrix(x%*%Usparse, ncol=K)
  if(K==1) xtranssparse <- matrix(xtranssparse, ncol=1)
  if(!is.null(xte)){
    #xtetranssparse <- matrix(xte%*%Usparse, ncol=K)
    xtetranssparse <- xte%*%Usparse
    if(K==1) xtetranssparse <- matrix(xtetranssparse,ncol=1)
    ypredsparsemat <- matrix(NA, ncol=K, nrow=nrow(xte))
    for(k in 1:K){
      ypredsparsemat[,k] <- Classify(matrix(xtranssparse[,1:k],ncol=k),matrix(xtetranssparse[,1:k],ncol=k),yclass)
    }
    obj <- (list(ypred=ypredsparsemat,discrim=Usparse,xproj=xtranssparse,xteproj=xtetranssparse, K=K, crits=penpca$crits,type=type,lambda=lambda,lambda2=lambda2,wcsd.x=wcsd.x,x=xorig,y=yclass))
    class(obj) <- "penlda"
    return(obj)
  }
  obj <- (list(discrim=Usparse,xproj=xtranssparse,K=K, crits=penpca$crits,type=type,lambda=lambda,lambda2=lambda2,wcsd.x=wcsd.x,x=xorig,y=yclass))
  class(obj) <- "penlda"
  return(obj)
}

predict.penlda <- function(object,xte,...){
  # First need to standardize the training and test data sets as needed.
  meanvec <- apply(object$x,2,mean)
  if(is.null(object$wcsd.x)){
    xte <- scale(xte,center=meanvec,scale=FALSE)
    x <- scale(object$x,center=meanvec,scale=FALSE)
  }
  if(!is.null(object$wcsd.x)){
    xte <- scale(xte,center=meanvec,scale=object$wcsd.x)
    x <- scale(object$x, center=meanvec,scale=object$wcsd.x)
  }
  # Done standardizing.
  # Now perform classification.
#  if(object$K==1){
#    ypred <- Classify(matrix(x%*%object$discrim,ncol=1), matrix(xte%*%object$discrim,ncol=1), object$y)
#    return(list(ypred=ypred))
#  }
  ypredsparsemat <- matrix(NA, nrow=nrow(xte), ncol=object$K)
  for(k in 1:object$K){
    ypredsparsemat[,k] <- Classify(matrix(x%*%object$discrim[,1:k],ncol=k), matrix(xte%*%object$discrim[,1:k],ncol=k), object$y)
  }
  return(list(ypred=ypredsparsemat))
}

plot.penlda <- function(x,...){
  K <- x$K
  par(mfrow=c(1,K))
  for(k in 1:K) plot(x$discrim[,k], main=paste("Discriminant ", k, sep=""),xlab="Feature Index", ylab="")
}

print.penlda <- function(x,...){
  cat("Number of discriminant vectors: ", x$K, fill=TRUE)
  K <- x$K
  for(k in 1:K){
    cat("Number of nonzero features in discriminant vector ", k, ":", sum(x$discrim[,k]!=0),fill=TRUE)
  }
  if(K>1) cat("Total number of nonzero features: ", sum(apply(x$discrim!=0, 1, sum)!=0),fill=TRUE)
  cat(fill=TRUE)
  cat("Details:", fill=TRUE)
  cat("Type: ", x$type, fill=TRUE)
  if(x$type=="standard") cat("Lambda: ", x$lambda, fill=TRUE)
  if(x$type=="ordered"){
    cat("Lambda: ", x$lambda,fill=TRUE)
    cat("Lambda2: ", x$lambda2, fill=TRUE)
  }
}
