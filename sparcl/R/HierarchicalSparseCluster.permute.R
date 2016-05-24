`HierarchicalSparseCluster.permute` <-
function(x,  nperms=10, wbounds=NULL, dissimilarity=c("squared.distance","absolute.value"), standardize.arrays=FALSE){
  dissimilarity <- match.arg(dissimilarity)
  if(is.null(wbounds)) wbounds <- seq(1.1, sqrt(ncol(x))*.7, len=10)
  if(min(wbounds)<=1) stop("Cannot have wbounds <= 1")
  if(length(wbounds)<2) stop("Wbounds should be a vector with at least 2 elements.")
  tots <- rep(NA,length(wbounds))
  permtots <- matrix(NA, nrow=length(wbounds), ncol=nperms)
  nnonzerows <- rep(NA,length(wbounds))
  if(standardize.arrays){
    #x <- t(scale(t(x),T,T))
    x <- sweep(x,1,apply(x,1,mean,na.rm=TRUE),"-")
    x <- sweep(x,1,apply(x,1,sd,na.rm=TRUE),"/")   
  }
  cat("Running sparse hierarchical clustering on unpermuted data",fill=TRUE)
  for(i in 1:length(wbounds)){ 
    cat(i,fill=FALSE)
    if(i==1) out <- HierarchicalSparseCluster(x,  wbound=wbounds[i],silent=TRUE,dissimilarity=dissimilarity)
    if(i>1) out <- HierarchicalSparseCluster(x=NULL,dists=out$dists,wbound=wbounds[i], silent=TRUE,dissimilarity=dissimilarity)
    nnonzerows[i] <- sum(out$ws!=0)
    tots[i] <- out$crit
  }
  cat(fill=TRUE)
  cat("Running sparse hierarchical clustering on permuted data",fill=TRUE)
  permdists <- out$dists
  for(k in 1:nperms){
    cat("Permutation ", k, " of ", nperms,fill=TRUE)
    # Oooohhhh.. It turns out that rather than permuting the columns of x and then computing a dist matrix, we can simply permute
    #  the columns of the (n choose 2)xp dist matrix.
    for(j in 1:ncol(permdists)) permdists[,j] <- sample(permdists[,j])
    for(i in 1:length(wbounds)){
      cat(i,fill=FALSE)
      perm.out <- HierarchicalSparseCluster(x=NULL, dists=permdists,wbound=wbounds[i], silent=TRUE,dissimilarity=dissimilarity)
      permtots[i,k] <- max(perm.out$crit)
    }
    cat(fill=TRUE)
  }
  gaps <- (log(tots)-apply(log(permtots),1,mean))
  out <- list(tots=tots, permtots=permtots, nnonzerows=nnonzerows, gaps=gaps, sdgaps=apply(log(permtots),1,sd), wbounds=wbounds, bestw=wbounds[which.max(gaps)], dists=out$dists)
  class(out) <- "HierarchicalSparseCluster.permute"
  cat(fill=TRUE)
  return(out)
}

print.HierarchicalSparseCluster.permute <- function(x,...){
  cat("Tuning parameter selection results for Sparse Hierarchical Clustering:", fill=TRUE)
  mat <- round(cbind(x$wbounds, x$nnonzerows, x$gaps, x$sdgaps),4)
  dimnames(mat) <- list(1:length(x$wbounds), c("Wbound", "# Non-Zero W's", "Gap Statistic", "Standard Deviation"))
  print(mat, quote=FALSE)
  cat("Tuning parameter that leads to largest Gap statistic: ", x$bestw, fill=TRUE)
}

plot.HierarchicalSparseCluster.permute <- function(x,...){
  plot(x$nnonzerows, x$gaps, log="x", main="Gap Statistics", xlab="# Non-zero Wj's", ylab="")
  lines(x$nnonzerows, x$gaps)
}
