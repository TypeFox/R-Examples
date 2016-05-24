`HierarchicalSparseCluster` <-
function(x=NULL, dists=NULL, method=c("average", "complete", "single","centroid"), wbound=NULL,niter=15,dissimilarity=c("squared.distance", "absolute.value"), uorth=NULL, silent=FALSE, cluster.features=FALSE, method.features=c("average", "complete", "single","centroid"),output.cluster.files=FALSE, outputfile.prefix="output",genenames=NULL, genedesc=NULL,standardize.arrays=FALSE){
  method <- match.arg(method)
  method.features <- match.arg(method.features)
  dissimilarity <- match.arg(dissimilarity)
  if (is.null(x) && is.null(dists)) stop("x or dists must be given!!!")
  xorig <- x
  if(standardize.arrays){
    if(is.null(x)) stop("Cannot standardize arrays if x not given.")
    dists <- NULL
    x <- sweep(x,1,apply(x,1,mean,na.rm=TRUE),"-")
    x <- sweep(x,1,apply(x,1,sd,na.rm=TRUE),"/")
  }
  if(is.null(dists)){
    xnona <- x
    xnona[is.na(x)] <- 0
    dists <- matrix(distfun(xnona), ncol=ncol(x)) # Rob's Fortran code!!! - no missing values please
    if(sum(is.na(x))>0){
      xbinary <- matrix(1, nrow=nrow(x),ncol=ncol(x))
      xbinary[is.na(x)] <- 0
      mult <- matrix(multfun(xbinary),ncol=ncol(x)) # mult equals 1 if neither elt is NA, and 0 if one or both is NA
      if(dissimilarity=="squared.distance"){
        dists <- sweep(dists,1,sqrt(ncol(dists)/apply(mult!=0,1,sum)),"*")
      } else if (dissimilarity=="absolute.value"){
        dists <- sweep(dists,1,ncol(dists)/apply(mult!=0,1,sum),"*")
      }
      dists[mult==0] <- 0
      mult <- NULL
      xbinary <- NULL
      xnona <- NULL
    }
  }
  if(is.null(wbound)) wbound <- .5*sqrt(ncol(dists))
  if(wbound<=1) stop("Cannot have wbound <= 1")
  if (dissimilarity == "squared.distance") out <- GetUW(dists^2, wbound, niter = niter, uorth = uorth, silent = silent)
  if (dissimilarity == "absolute.value")  out <- GetUW(dists, wbound, niter = niter, uorth = uorth, silent = silent)
  out <- list(hc = hclust(as.dist(out$u), method = method), ws = out$w, u = out$u, crit = out$crit, dists = dists, uorth = uorth, wbound = wbound)

  if(cluster.features){
    if(is.null(x)) stop("Cannot cluster features unless x is given.")
    rho=cor(xorig[,out$ws!=0],use="pairwise.complete.obs")
    out2 <- hclust(as.dist(2*(1-rho)), method=method.features)
    out$hc.features <- out2
  }

  if(output.cluster.files){
    if(is.null(x)) stop("Cannot output files unless x is given.")
    output.cluster.files.fun(xorig,out,outputfile.prefix,genenames=genenames,genedesc=genedesc)
  }
  if(!silent) cat(fill=TRUE)
  class(out) <- "HierarchicalSparseCluster"
  return(out)
}


plot.HierarchicalSparseCluster <- function(x,...){
  par(mfrow=c(1,2))
  plot(x$hc,xlab="",ylab="",sub="", main="Sparse Clustering", labels=rep("", nrow(x$u)))
  plot(x$ws, main=paste("Wbound is ", sep="", round(x$wbound,3)), xlab="Feature Index", ylab="Wj")
}



print.HierarchicalSparseCluster <- function(x,...){
  cat("Wbound is ", x$wbound, ":", fill=TRUE)
  cat("Number of non-zero weights: ", sum(x$ws!=0), fill=TRUE)
  cat("Sum of weights: ", sum(x$ws), fill=TRUE)
  cat(fill=TRUE)
}
