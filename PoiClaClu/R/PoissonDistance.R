PoissonDistance <-
function(x,beta=1,type=c("mle","deseq","quantile"),transform=TRUE, alpha=NULL,perfeature=FALSE){
  type <- match.arg(type)
  if(!transform && !is.null(alpha)) stop("You have asked for NO transformation but have entered alpha.")
  if(transform && !is.null(alpha)){
    if(alpha>0 && alpha <= 1) x <- x^alpha
    if(alpha<=0 || alpha>1) stop("alpha must be between 0 and 1")
  }  
  if(transform && is.null(alpha)){
    alpha <- FindBestTransform(x)
    x <- x^alpha
  }
  dd <- matrix(0, nrow=nrow(x), ncol=nrow(x))
  ddd <- NULL
  if(perfeature) ddd <- array(0, dim=c(nrow(x), nrow(x), ncol(x)))
  for(i in 2:nrow(dd)){
    xi <- x[i,]
    for(j in 1:(i-1)){
      xj <- x[j,]
      n <- NullModel(x[c(i,j),],type=type)$n
      ni <- n[1,]
      nj <- n[2,]
      di <- (xi+beta)/(ni+beta)
      dj <- (xj+beta)/(nj+beta)
      dd[i,j] <- sum(ni+nj-ni*di-nj*dj+xi*log(di)+xj*log(dj))
      if(perfeature) ddd[j,i,] <- ddd[i,j,] <- ni+nj-ni*di-nj*dj+xi*log(di)+xj*log(dj)
    }
  }
  save <- list(dd=as.dist(dd+t(dd)), alpha=alpha, x=x, ddd=ddd, alpha=alpha, type=type)
  class(save) <- "poidist"
  return(save)
}

print.poidist <- function(x,...){
  if(!is.null(x$alpha)) cat("Value of alpha used to transform data: ", x$alpha, fill=TRUE)
  if(is.null(x$alpha)) cat("No transformation performed.",fill=TRUE)
  cat("This type of normalization was performed:", x$type, fill=TRUE)
  cat("Dissimilarity computed for ", nrow(x$x), " observations.", fill=TRUE)
}

