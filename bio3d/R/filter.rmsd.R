filter.rmsd <- function(xyz=NULL, rmsd.mat=NULL, cutoff=0.5, fit=TRUE, verbose=TRUE,
                        inds=NULL, ...) {

  # k<-filter.rmsd(xyz=pdbs$xyz, cutoff=0.5)
  # k<-filter.rmsd(rmsd.mat=k$rmsd.mat, cutoff=2.0)
  
  if(is.null(rmsd.mat)) {
    if(is.null(xyz))
      stop("Must provide either a 'xyz' matrix or RMSD matrix 'rmsd.mat'")
    if(is.list(xyz))
      xyz=xyz$xyz
    if(is.null(inds)) {
       gaps <- gap.inspect(xyz)
       inds <- gaps$f.inds
    }
    rmsd.mat <- rmsd( xyz, a.inds=inds, fit=fit, ... )
  }
  
  r.d  <- as.dist(rmsd.mat)
  tree <- hclust(r.d)

  h <- cutoff
  n <- nrow(tree$merge) + 1
  k <- integer(length(h))
  k <- n + 1 - apply(outer(c(tree$height, Inf), h, ">"),2, which.max)

  if(verbose)
    cat("filter.rmsd(): N clusters @ cutoff = ", k, "\n")
  
  #ans <- as.vector(.Call("R_cutree", tree$merge, k, PACKAGE = "stats"))
  ans <- as.vector(cutree(tree, k))
  
  cluster.rep <- NULL
  for(i in 1:k) {
    ind <- which(ans==i)
    if (length(ind) == 1) {
      cluster.rep <- c(cluster.rep, ind)
    } else {
      cluster.rep <- c(cluster.rep,
                     ind[ which.min( colSums(rmsd.mat[ind,ind]) ) ])
    }
  }

  return(list(ind=cluster.rep, tree=tree, rmsd.mat=rmsd.mat))
}

