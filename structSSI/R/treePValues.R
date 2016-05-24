treePValues <- function(tree, abundances, environments){

  igraphTree <- graph.edgelist(tree)

  treePValues <- vector(length = length(V(igraphTree)))
  names(treePValues) <- V(igraphTree)$name

  for(vertex in V(igraphTree)){

    ## first, aggregate data descending ##
    ## from that node ##
    graphDiameter <- diameter(igraphTree)
    curVertexName <- V(igraphTree)[vertex]$name
    descendants <- neighborhood(igraphTree, nodes = curVertexName, mode = "out", order = graphDiameter)
    subtree <- induced.subgraph(igraphTree, vids = descendants[[1]])

    ## some nodes don't have descendants that are in the OTU table,
    ## so we can't consider them in the p-value calculations.
    allDescendantsNames <- V(subtree)$name
    namesInTipsIndex <- which(allDescendantsNames %in% rownames(abundances))

    if(length(namesInTipsIndex) > 0){

      subOtuTable <- abundances[V(subtree)$name[namesInTipsIndex], , drop=F]
      aggregateData <- colSums(subOtuTable)
      aggregateDataWithLabel <- data.frame(abund = as.vector(aggregateData),
                                           type = environments)

      ## fit the LM to obtain p-values ##

      abundanceModel <- summary(lm(abund ~ type, data = aggregateDataWithLabel))
      treePValues[curVertexName] <- pf(abundanceModel$fstatistic[1],
                                       abundanceModel$fstatistic[2],
                                       abundanceModel$fstatistic[3],
                                       lower.tail = FALSE)
    }
  }
  return(treePValues)
}
