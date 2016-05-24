Classify.cv <-
function(x,y,rhos=NULL,beta=1,nfolds=5,type=c("mle","deseq","quantile"),folds=NULL,transform=TRUE, alpha=NULL, prior=NULL){
  type <- match.arg(type)
  if(!transform && !is.null(alpha)) stop("You have asked for NO transformation but have entered alpha.")
  if(transform && is.null(alpha)) alpha <- FindBestTransform(x)
   if(transform){
     if(alpha<=0 || alpha>1) stop("alpha must be between 0 and 1")
     x <- x^alpha
  }
  if(is.null(rhos)){
    ns <- NullModel(x,type=type)$n
    uniq <- sort(unique(y))
    maxrho <- rep(NA, length(uniq))
    for(k in 1:length(uniq)){
      a <- colSums(x[y==uniq[k],])+beta
      b <- colSums(ns[y==uniq[k],])+beta
      maxrho[k] <- max(abs(a/b-1)*sqrt(b),na.rm=TRUE)
    }
    rhos <- seq(0, max(maxrho,na.rm=TRUE)*(2/3), len=30)
  }
  if(is.null(folds)) folds <- balanced.folds(y,nfolds=nfolds)
  nfolds <- length(folds)
  errs <- nnonzero <- matrix(NA, nrow=nfolds, ncol=length(rhos))
  for(i in 1:nfolds){
    cat(i,fill=FALSE)
    tr <- -folds[[i]]
    te <- folds[[i]]
    out <- Classify(x[tr,],y[tr],x[te,],rhos=rhos,beta=beta,type="quantile", prior=prior, transform=FALSE) # Have already power-transformed x, so don't need to do it again!!!
    for(j in 1:length(rhos)){      
      errs[i,j] <- sum(out[[j]]$ytehat!=y[te])
      nnonzero[i,j] <- sum(colSums(out[[j]]$ds!=1)!=0)
    }
  }
  cat(fill=TRUE)
  save <- list(errs=errs, bestrho=rhos[max(which(colMeans(errs)==min(colMeans(errs))))], rhos=rhos, nnonzero=nnonzero,folds=folds, alpha=alpha,type=type)
  class(save) <- "poiclacv"
  return(save)
}

print.poiclacv <- function(x,...){
  if(!is.null(x$alpha)) cat("Value of alpha used to transform data: ", x$alpha, fill=TRUE)
  if(is.null(x$alpha)) cat("No transformation performed.",fill=TRUE)
  cat("Rho values considered: ", round(x$rhos,3), fill=TRUE)
  cat("Number of CV folds performed: ", length(x$folds), fill=TRUE)
  cat("Type of normalization performed: ", x$type, fill=TRUE)
  nnonzero.mean <- round(apply(x$nnonzero, 2, mean),3)
  err.mean <- round(apply(x$errs, 2, mean),3)
  cat(fill=TRUE)
  cat("CV results:", fill=TRUE)
  mat <- cbind(round(x$rhos,3), err.mean,nnonzero.mean)
  mat <- rbind(c("Rho", "Errors", "Num. Nonzero Features"), mat)
  write.table(mat, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}



