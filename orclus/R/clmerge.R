### sub functions to merge clusters during iterations
# compute projected new subspace and projected energy merge of two clusters
compute.projen <- function(ij, x, act.clustering){
  if (is.null(ij[[2]])){
    pair        <- ij[[1]]
    mergedx     <- x[c(which(act.clustering$cluster == pair[1]), which(act.clustering$cluster == pair[2])),] # merged cluster
    if(is.vector(act.clustering$subspaces[[1]])) actu.l <- 1 # actual subspace dimension
    if(is.matrix(act.clustering$subspaces[[1]])) actu.l <- dim(act.clustering$subspaces[[1]])[2] 
    subspace    <- eigen(cov(mergedx))$vectors   # subspace of merged cluster
    subspace    <- subspace[,(dim(subspace)[2]+1-actu.l):dim(subspace)[2]]
    cen         <- apply(mergedx %*% subspace, 2, mean) # centers of merged cluster in its subspace
    centpro     <- t(apply(mergedx %*% subspace, 1, function(z, zen) return(z-zen), zen = cen)) # merged cluster in subspace shifted by cluster mean 
    projen      <- sum(apply(centpro, 1, function(z) return(sum(z^2)))) / nrow(centpro)  # projected energy: average distances from centers in subspace
    }
  if (!is.null(ij[[2]])){
    pair        <- ij[[1]]
    subspace    <- ij[[2]]
    projen      <- ij[[3]]
    }    
  result        <- list(pair, subspace, projen)
  names(result) <- c("pair", "subspace", "projen")  
  return(result)
  }

# update clustering result after merging 2 clusters  
update.actclu <- function(win.pair, newsubsp, act.clustering, x){
   # assign new common subspace of both clusters
   act.clustering$subspaces[[win.pair[1]]] <- newsubsp
   # merge cluster labels 
   act.clustering$cluster[act.clustering$cluster == win.pair[2]] <- win.pair[1]
   act.clustering$cluster[act.clustering$cluster > win.pair[2]] <- act.clustering$cluster[act.clustering$cluster > win.pair[2]] -1
   # compute new cluster sizes
   act.clustering$size <- table(act.clustering$cluster) 
   # compute new cluster means in original space
   newmeans <- by(x, act.clustering$cluster, colMeans)
   newcenters <- NULL
   for (i in 1:length(table(act.clustering$cluster))) newcenters <- rbind(newcenters, newmeans[[i]])
   act.clustering$centers <- newcenters
   # reorder the list of clusterspecific subspaces
   for(i in win.pair[2]:length(act.clustering$subspaces)){
      if(i > win.pair[2]) act.clustering$subspaces[[i-1]] <- act.clustering$subspaces[[i]]
      }
   length(act.clustering$subspaces) <- length(act.clustering$subspaces)-1
   return(act.clustering)   
  }  

# update the inter-cluster results to be re-used for next merge step...
update.mergelist <- function(win.pair, mergelist, act.clustering){
  i <- win.pair[1]
  j <- win.pair[2]  
  #recode clusterids in mergelist
  change.pairs <- function(mlelem, i = win.pair[1], j = win.pair[2]){
    id1 <- mlelem[[1]][1]
    id2 <- mlelem[[1]][2]
    if(id1 == j | id2 == j) {mlelem[[1]] <- c(0,0); id1 <- 0; id2 <- 0}
    if(id1 > j) mlelem[[1]][1] <- id1 - 1    
    if(id2 > j) mlelem[[1]][2] <- id2 - 1    
    return(mlelem)
    }
  mergelist.old <- lapply(mergelist, change.pairs, i = win.pair[1], j = win.pair[2])
  # generate new mergelist and copy old calculations that are still valid
  ml.new <- list()
  i <- 1
  for(j in 1:length(mergelist.old)){
  if (sum(mergelist.old[[j]][[1]] == 0) != 2){ 
    ml.new[[i]] <- mergelist.old[[j]] 
    i <- i+1
    }
  }
  # remove unusual results of the merged new cluster  
  for(i in 1:length(ml.new)){  
    if(ml.new[[i]][[1]][1] == win.pair[1] | ml.new[[i]][[1]][2] == win.pair[1]){
      ml.new[[i]][[2]] <- ml.new[[i]][[3]] <- NULL
      length(ml.new[[i]]) <- 3
      names(ml.new[[i]]) <- c("pair", "subspace", "projen")
      }
     }  
  return(ml.new)  
  }

# iteratively merge two clusters until new number of clusters is reached
clmerge <- function(x, act.clustering, knew){
  # only merge clusters if actual number is not reached
  if (length(act.clustering$size) > knew){
    # potential cluster pairs to be merged
    clpairs   <- t(combn(length(table(act.clustering$cluster)), 2))
    #initialize list with potential merging results
    mergelist <- list()
    length(mergelist) <- nrow(clpairs)
    for (i in 1:length(mergelist)){
      mergelist[[i]] <- list()
      length(mergelist[[i]]) <- 3
      mergelist[[i]][[1]]   <- clpairs[i,]
      names(mergelist[[i]]) <- c("pair", "subspace", "projen")
      }
    for (mergesteps in 1:(length(act.clustering$size) - knew)){
      if (length(table(act.clustering$cluster)) > knew){
         # compute subspaces and projected energies for any potential merges
         mergelist <- lapply(mergelist, compute.projen, x = x, act.clustering = act.clustering)
         # find cluster pair with minimal projected energy     
         projens <- which.min(sapply(mergelist, function(z) return(z[[3]])))   
         if (length(projens) > 1) projens <- sample(projens, 1) # random choice if more than one winner
         win.pair <- mergelist[[projens]]$pair
         newsubsp <- mergelist[[projens]]$subspace
         # update the clustering result for the merged clusters    
         act.clustering <- update.actclu(win.pair, newsubsp, act.clustering, x)
         # update list of results for the next merge step 
         mergelist <- update.mergelist(win.pair, mergelist, act.clustering) 
         }
      }   
    }        
  return(act.clustering)  
  }

