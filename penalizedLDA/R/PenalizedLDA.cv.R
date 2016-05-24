PenalizedLDA.cv <-
function(x, y, lambdas=NULL, K=NULL, nfold=6, folds=NULL, type="standard", chrom=NULL, lambda2=NULL){
  if(is.null(folds)){
    folds <- balanced.folds(y, nfold)
  } else {
    nfold <- length(folds)
  }
  if(sum(1:length(unique(y)) != sort(unique(y)))>0) stop("y must be a numeric vector, with values as follows: 1, 2, ....")
  if(sum(is.na(x))>0 || sum(is.na(y))>0) stop("No missing values allowed!!!")
  

  if(is.null(lambdas)) lambdas <- seq(0, 5, len=10)
  if(length(unique(y))==2 || !is.null(K)){
    if(is.null(K)) K <- 1
    err <- matrix(NA, nrow=nfold, ncol=length(lambdas))
    nzerobetas <- matrix(NA, nrow=nfold, ncol=length(lambdas))
    for(fold in 1:nfold){
      cat("Fold ", fold, fill=TRUE)
      xtr <- x[-folds[[fold]],]
      ytr <- y[-folds[[fold]]]
      xte <- x[folds[[fold]],]
      yte <- y[folds[[fold]]]
      # Doing some centering and computations now so that don't need to do it for each value of lambda #
      if(length(ytr)<200){
        wcsd.x <- wcsd.matrix(xtr,MakeYMat(ytr))#apply(xtr, 2, wcsd, ytr)
      } else {
        wcsd.x <- apply(xtr,2,wcsd,ytr)
      }
      if(min(wcsd.x)==0) stop("Some features have within-class standard deviation equal to zero.")
#      wcsd.x <- wcsd.x + quantile(wcsd.x, .05)
      xte <- scale(xte, center=apply(xtr,2,mean), scale=wcsd.x)
      xtr <- scale(xtr, T, scale=wcsd.x)
      ymat <- MakeYMat(ytr)
      # Done these computations #
      for(i in 1:length(lambdas)){
        cat(i,fill=FALSE)
        out <- PenalizedLDA(x=xtr, y=ytr, xte=xte, lambda=lambdas[i], K=K, type=type, chrom=chrom, lambda2=lambda2, standardized=TRUE, ymat=ymat)
        if(length(out$ypred)==0) out$ypred <- rep(-1, length(yte))
        if(K==1) err[fold,i] <- sum(out$ypred!=yte)
        if(K>1) err[fold,i] <- sum(out$ypred[,K]!=yte)
        nzerobetas[fold,i] <- sum(apply(abs(out$discrim),1,sum)!=0)
      }
    }
    errmean <- apply(err, 2, mean)
    errse <- apply(err,2,sd)/sqrt(nfold)
    bestlambda.1se <- max(lambdas[errmean <= (min(errmean)+errse)])
    obj <- list(errs=errmean, lambdas=lambdas, bestlambda=mean(lambdas[which(errmean==min(errmean))]), 
             nnonzero=apply(nzerobetas, 2, mean), bestK=K,folds=folds,bestlambda.1se=bestlambda.1se)
    class(obj) <- "penldacv"
    return(obj)
  }
  # Now, we're in the case where there are more than 2 classes, and K hasn't been specified
  # So we'll cross-validate over K
  Ks <- 1:(length(unique(y))-1)
  err <- array(NA, dim=c(nfold, length(lambdas), length(Ks)))
  nzerobetas <- array(NA, dim=c(nfold, length(lambdas), length(Ks)))
  for(fold in 1:nfold){
    cat("Fold ", fold, fill=TRUE)
    xtr <- x[-folds[[fold]],]
    ytr <- y[-folds[[fold]]]
    xte <- x[folds[[fold]],]
    yte <- y[folds[[fold]]]
    for(i in 1:length(lambdas)){
      cat(i,fill=FALSE)
      out <- PenalizedLDA(x=xtr, y=ytr, xte=xte, lambda=lambdas[i], K=max(Ks), type=type, chrom=chrom, lambda2=lambda2)
      for(k in 1:length(Ks)){
        xtranssparse <- matrix(x%*%out$discrim[,1:Ks[k]], ncol=Ks[k])
        xtetranssparse <- matrix(xte%*%out$discrim[,1:Ks[k]], ncol=Ks[k])
        ypredsparse <- Classify(xtranssparse,xtetranssparse,y)
        if(length(ypredsparse)==0) ypredsparse <- rep(-1, nrow(xte))
        err[fold,i,k] <- sum(ypredsparse!=yte)
        if(Ks[k]>1) nzerobetas[fold,i,k] <- sum(apply(abs(out$discrim[,1:Ks[k]]), 1, sum)!=0)
        if(Ks[k]==1) nzerobetas[fold,i,k] <- sum(out$discrim[,1]!=0)
      }
    }
  }
  errmean <- apply(err, c(2,3), mean)
  errse <- apply(err, c(2,3), sd)/sqrt(nfold)
  nnonzerobeta <- apply(nzerobetas, c(2,3), mean)
  mat <- which(errmean==min(errmean), arr.ind=TRUE)
  bestK <- Ks[min(mat[,2])]
  whichers <- which(mat[,2]==min(mat[,2]))
  bestlambda <- max(lambdas[mat[whichers,1]])
  bestlambda.1se <- max(lambdas[errmean[,bestK] <= (min(errmean[,bestK])+errse[,bestK])])
  cat("Best K is ", bestK, fill=TRUE)
  obj <- list(errs=errmean, nnonzero=nnonzerobeta, Ks=Ks, lambdas=lambdas,
              bestK=bestK, bestlambda=bestlambda,folds=folds,bestlambda.1se=bestlambda.1se)
  class(obj) <- "penldacv"
  return(obj)
}


print.penldacv <- function(x,...){
  cat("Cross-validation Results:", fill=TRUE)
  cat("Values of Lambda considered: ", round(x$lambdas,3), fill=TRUE)
  if(!is.null(x$Ks)){
    cat("Values of K considered: ", x$Ks, fill=TRUE)
    mat.errs <- x$errs
    cat("Matrix of Cross-Validation Errors:", fill=TRUE)
    write.table(round(mat.errs,3),sep="\t",col.names=paste("K=",x$Ks,sep=""), row.names=paste("Lambda=",round(x$lambdas,3),sep=""), quote=FALSE)
    mat.nnonzero <- x$nnonzero
    cat("Matrix of Number of Nonzero Features:", fill=TRUE)
    write.table(round(mat.nnonzero,3),sep="\t",col.names=paste("K=",x$Ks,sep=""), row.names=paste("Lambda=",round(x$lambdas,3),sep=""), quote=FALSE)
    cat("Best Lambda: ", x$bestlambda, " Best K : ", x$bestK, fill=TRUE)
  } else if(is.null(x$Ks)){
    cat("Used only 1 value of K: ", x$bestK, fill=TRUE)
    cat("Mean CV Errors: ", round(x$errs,3), fill=TRUE)
    cat("Mean Number of Nonzero Features: ", round(x$nnonzero, 3), fill=TRUE)
    cat("Value of Lambda with lowest CV error: ", round(x$bestlambda,3), fill=TRUE)
  }
}

plot.penldacv <- function(x,...){
  if(is.null(x$Ks) || length(x$Ks)==1){
    plot(x$nnonzero, x$errs, main="CV Error", xlab="Num. Nonzero Features", ylab="CV Error")
  } else {
    for(k in x$Ks){
      if(k==1) plot(x$nnonzero[,k], x$errs[,k], main="CV Error", xlab="Num. Nonzero Features", ylab="CV Error", type="l", xlim=range(x$nnonzero),ylim=range(x$errs))
      if(k>1) points(x$nnonzero[,k], x$errs[,k], col=k, type="l")
      legend("topright", pch=15, col=(x$Ks), paste(x$Ks, " Disc. Vecs Used", sep=""))
    }
  }
}
