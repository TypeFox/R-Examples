filter.identity <- function(aln=NULL, ide=NULL, cutoff=0.6, verbose=TRUE, ...) {
                            
  
 #k<-filter.identity(aln,cutoff=0.4)
 #aln$id[k$ind]
 #k<-filter.identity(ide=k$ide,cutoff=0.6)
 #plot(k$tree, axes = FALSE, ylab="%identity")
 #axis(2,labels =c(1,0.8,0.6,0.4))
 #abline(h=0.6)

  if(is.null(ide)) {
    if(is.null(aln)) 
      stop("Must provide either an alignment 'aln' or identity matrix 'ide'")
    ide  <- seqidentity(aln, ...)
  }
  i.d  <- as.dist(1-ide)
  tree <- hclust(i.d)

  h <- 1 - cutoff
  n <- nrow(tree$merge) + 1
  k <- integer(length(h))
  k <- n + 1 - apply(outer(c(tree$height, Inf), h, ">"),2, which.max)
  if(verbose)
    cat("filter.identity(): N clusters @ cutoff = ", k, "\n")
  
  #ans <- as.vector(.Call("R_cutree", tree$merge, k, PACKAGE = "stats"))
  ans <- as.vector(cutree(tree, k))

  cluster.rep <- NULL
  for(i in 1:k) {
    ind <- which(ans==i)
    if (length(ind) == 1) {
      cluster.rep <- c(cluster.rep, ind)
    } else {
      cluster.rep <- c(cluster.rep,
                     ind[ which.max( colSums(ide[ind,ind]) ) ]) # max similarity
    }
  }
  return(list(ind=cluster.rep, tree=tree, ide=ide))
}

