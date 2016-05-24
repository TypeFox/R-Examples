getTipList <- function(phy) {
  tipList <- cbind(phy$edge, phy$edge[,2] %in% phy$edge[,1], rep(0, dim(phy$edge)[1]))
  tips <- which(tipList[,3] == 0)
  tipList[tips,3] <- "tip"
  tipList[tips,4] <- phy$tip.label[as.numeric(tipList[tips,2])]
  ints <- which(tipList[,3] == 1)
  tipList[ints,3] <- "internal"
  return(data.frame(tipList, stringsAsFactors=F))
}

whichEdge <- function(phy, taxa) {
  tipList <- getTipList(phy)
  nodes <- tipList[tipList[, 4] %in% taxa, 1]
  foundBranch <- FALSE
  while(foundBranch==FALSE) {
    for(i in as.numeric(unique(nodes))){
      leaves <- node.leaves(phy, i)
      if(all(leaves %in% taxa) && all(taxa %in% leaves)){
        foundBranch <- TRUE
        return(as.numeric(i))  
      }
    }
    nodes <- tipList[which(tipList[,2] %in% nodes), 1]
  }
}

node.leaves<-function(phy, node) {
  n<-length(phy$tip.label)
  if(node<= n) return(phy$tip.label[as.numeric(node)])
    l<-character();
  d<-node.offspring(phy, node);
  for(j in d) {
    if(j <= n) l<-c(l, phy$tip.label[as.numeric(j)])
    else l<-c(l, node.leaves(phy, j));
  }
  return(l);
}

node.offspring<-function(phy, node) {
  r<-which(phy$edge[,1]==node)
  return(phy$edge[r,2])
}


WhatToDoWithDuplicateEdgeNames <- function(edgeLabels, duplicateEdgeLabels){
  if(duplicateEdgeLabels == "recent")
    return(names(edgeLabels[length(edgeLabels)]))
  if(duplicateEdgeLabels == "oldest")
    return(names(edgeLabels[1]))
  if(duplicateEdgeLabels == "combined")
    return(paste(names(edgeLabels), sep="", collapse="."))
}

MakeEdgeLabels <- function(MyHiers, label="all", missingData=NULL, duplicateEdgeLabels="oldest"){
  MyHiers <- RemoveNAFiles(MyHiers)
  nodeList <- NodeLabelList(MyHiers, label="all", missingData=missingData)
  if(length(nodeList) == 0)
    stop("Node Labels can not be created, because hierarchy information doesn't overlap")
  phy <- MakeHierarchyTree(MyHiers, missingData=missingData, includeNodeLabels=FALSE)
  tipList <- getTipList(phy)
  edges <- c(lapply(nodeList, whichEdge, phy=phy), recursive=T)
  if(any(duplicated(names(edges))))
    edges <- edges[-which(duplicated(names(edges)))]  
  if(any(duplicated(edges))){
    for(i in unique(edges[which(duplicated(edges))])){
      duplicateEdges <- edges[which(edges == i)]
      names(edges)[which(edges == i)] <- WhatToDoWithDuplicateEdgeNames(duplicateEdges, duplicateEdgeLabels)
    }
  }
  for(i in sequence(length(edges))){
    tipList[which(tipList[,2] == edges[i]),4] <- names(edges)[i]
  }
  justInts <- which(tipList[,3] == "internal")   #to get row to use in ape
  names(justInts) <- tipList[tipList[,3] == "internal",4]  #associate taxon name with row
  if(any(names(justInts) == 0))
    justInts <- justInts[-which(names(justInts) == 0)]
  return(justInts)
}













