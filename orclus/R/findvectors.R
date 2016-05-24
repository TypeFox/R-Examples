### sub functions for clusterwise subspace calculation:
# subspace for a single cluster
compute.newsubs <- function(clid, x, act.clustering, dimen){
  clusx      <- x[which(act.clustering$cluster == clid),]
  evd        <- eigen(cov(clusx))
  if (sum(evd$values == 0) > dimen ) warning("Number of eigenvalues equal 0 is larger than the actual subspace dimension for at least one cluster. The choice of the subspace is arbitrary at this point.")
  subspace.k <- evd$vectors[, (ncol(clusx)+1-dimen) : ncol(clusx)]
  return(subspace.k)
  }

# compute subspaces for all clusters   
findvectors <- function(x, act.clustering, dimen){
  clids <- 1:length(table(act.clustering$cluster))
  subspaces <- lapply(clids, compute.newsubs, x = x, act.clustering = act.clustering, dimen = dimen)  
  return(subspaces)
  }     

        