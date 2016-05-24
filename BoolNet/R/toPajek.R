# Export a state table in <stateGraph> to a Pajek graph
toPajek <- function (stateGraph, file="boolean.net", includeLabels=FALSE, ...) 
{
  args <- list(...)
  
  if (!is.null(args$attractorInfo))
  {
    warning("The parameter \"attractorInfo\" is deprecated. Use \"stateGraph\" instead!")
    stateGraph <- args$attractorInfo
  }
  
  if (!inherits(stateGraph,"TransitionTable"))
    stateGraph <- getTransitionTable(stateGraph)
  
  geneCols <- setdiff(colnames(stateGraph),c("attractorAssignment","transitionsToAttractor"))
  numGenes <- (length(geneCols)) / 2
  
  from <- apply(stateGraph[,1:numGenes,drop=FALSE],1,paste,collapse="")
  to <- apply(stateGraph[,((numGenes+1):(2*numGenes)),drop=FALSE],1,paste,collapse="")
  vertices <- unique(c(from,to))
  vertexIndices <- seq_along(vertices)
  names(vertexIndices) <- vertices
  
  from <- vertexIndices[from]
  to <- vertexIndices[to]

 
  sink(file)
  cat("*Vertices ", length(vertices), "\r\n", sep = "")
  if (includeLabels)
  {
    for (i  in seq_along(vertices))
      cat(i," \"",vertices[i],"\"\r\n",sep="")
  }
  cat("*Arcs\r\n")

  for (i in seq_along(from))
    cat(from[i]," ",to[i]," 1\r\n",sep="")
  sink()
}
