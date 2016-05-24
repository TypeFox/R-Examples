
# Pieces that are common to block models.

SBM.ID.rotation <- function(ID.labels, label.count=max(ID.labels)) {
  #ID.labels=c(3,4,4,1,1,2,5); label.count=max(ID.labels)
  
  #first pin down the keepers.
  assigned <- replaced <- rep(NA, label.count)
  replacements <- rep(NA, length(ID.labels))
  
  for (ii in 1:label.count) {
    if (is.na(replacements[ii])) {
      assigned[ID.labels[ii]] <- ii
      replaced[ii] <- ID.labels[ii]
      replacements[ID.labels==ID.labels[ii]] <- ii
    }
  }
  while (sum(is.na(assigned))>0) {
    #node <- min(which(is.na(replacements)))
    vacancy <- min(which(is.na(assigned)))
    newlabel <- min(which(is.na(replaced)))
    
    assigned[vacancy] <- newlabel
    replaced[newlabel] <- vacancy
    
    replacements[ID.labels==vacancy] <- newlabel
  }
  
  return (assigned)
}

SBM.rotate.bvector <- function(sbm.b.vector, rotation) {
  bvt <- cbind(make.edge.list.selfies(length(rotation)), sbm.b.vector)
  bvt[,1] <- rotation[bvt[,1]]; bvt[,2] <- rotation[bvt[,2]];
  bvt[,1:2] <- t(apply(rbind(bvt[,1:2]), 1, sort))
  bvt <- rbind(bvt[order(bvt[,1], bvt[,2]),])
  bvt[,3]
}

SBM.rotate.block <- function (sbm.block, rotation) {
  sbm.block[rotation, rotation]
}

MMSBM.ID.rotation <- function(node.prob.table, label.count=max(ID.labels)) {
  whichmax.sampler <- function(x) {
    items <- which(max(x)==x)
    ifelse (length(items)==1, items, sample(items,1))
  }
  ID.labels <- apply(node.prob.table, 2, whichmax.sampler)
  return(SBM.ID.rotation (ID.labels, label.count=max(ID.labels)))
}



