### sub functions to re-assign objects to clusters given the subspaces
# compute distance of an observation to a single cluster
compute.dist <- function(clid, obs, act.clustering){
  pro.obs    <- t(as.matrix(obs)) %*% act.clustering$subspaces[[clid]]  # projection of the observation
  pro.cen    <- act.clustering$centers[clid,] %*% act.clustering$subspaces[[clid]]  # projection of cluster centers
  distancek  <- sum((pro.obs-pro.cen)^2)
  return(distancek)
  }

# compute distances of a single observation to all clusters and return closest cluster
compute.newcluster <- function(obs, act.clustering){
  clids         <- 1:length(act.clustering$subspaces)
  dists         <- sapply(clids, compute.dist, obs = obs, act.clustering = act.clustering)
  newcluster    <- which.min(dists)
  if (length(newcluster) > 1) newcluster <- sample(newcluster, 1)
  return(newcluster)  
  }

# reassign all clusters and compute update elements for act.clustering object (cf. subspaces)
reassign <- function(x, act.clustering){
  recomputation <- FALSE 
  res <- list()
  res$subspaces <- act.clustering$subspaces
  
  # reassignment 
  res$cluster <- apply(x, 1, compute.newcluster, act.clustering = act.clustering)

  # compute new cluster centers
  res$means <- by(x, res$cluster,colMeans)

  # nonempty cluster indices 
  res$centers <- NULL
  for (i in 1:length(res$means)) res$centers <- rbind(res$centers, as.numeric(res$means[[i]]))
  res$size <- table(res$cluster)

  if (length(res$size) == 1) warning("All objects are in the same cluster!")

  # identify and remove clusters (centers and subspaces) with 0 or only 1 element 
  remove.subspaces  <- (1:length(act.clustering$subspaces))[-as.numeric(names(which(res$size > 1)))]
  keep.centers      <- res$size > 1
  if (length(remove.subspaces) > 0){

  warning("At least one empty or single element cluster removed during iteration process...")
  recomputation <- TRUE
  act.clustering$centers   <- res$centers[keep.centers,]
  for(i in remove.subspaces[length(remove.subspaces):1]) act.clustering$subspaces[[i]] <- NULL 
  }
  if (recomputation) res <- reassign(x, act.clustering)  
  
  return(res)
  }                                                                       
