# Plots a graph that visualizes the state transitions and attractor basins. <attractorInfo> is an object
# of class AttractorInfo. This requires the igraph package.
# If <highlightAttractors> is set, attractor edges are drawn bold.
# If <colorBasins> is true, each basin is drawn in a different color. 
# Colors can be provided in <colorSet>.
# <layout> specifies the graph layouting function.
# If <piecewise> is true, subgraphs are layouted separately.
# <basin.lty> and <attractor.lty> specify the line types used to draw states in the basins
# and in the attractors (if <highlightAttractor> is set).
# If <plotIt> is not set, only the igraph object is returned, but no graph is plotted.
# ... provides further graphical parameters for the plot.
# Returns an object of class igraph
plotStateGraph <- function(stateGraph,highlightAttractors=TRUE,colorBasins=TRUE,colorSet,
             drawLegend=TRUE,drawLabels=FALSE,layout=layout.fruchterman.reingold,
             piecewise=FALSE,basin.lty=2,attractor.lty=1,plotIt=TRUE,...)
{
  stopifnot(inherits(stateGraph,"AttractorInfo") || inherits(stateGraph,"TransitionTable") || 
            inherits(stateGraph,"SymbolicSimulation"))

  args <- list(...)
  
  if (!is.null(args$attractorInfo))
  {
    warning("The parameter \"attractorInfo\" is deprecated. Use \"stateGraph\" instead!")
    stateGraph <- args$attractorInfo
  }

  if (installed.packages()["igraph","Version"] < package_version("0.6"))
    bias <- 1
  else
    bias <- 0
  
  symbolic <- FALSE
  if (inherits(stateGraph,"AttractorInfo"))
  {
    stateGraph <- getTransitionTable(stateGraph)
  }
  else
  if (inherits(stateGraph,"SymbolicSimulation"))
  {
    symbolic <- TRUE
    
     if (is.null(stateGraph$graph))
      stop(paste("This SymbolicSimulation structure does not contain transition table information.",
                 "Please re-run simulateSymbolicModel() with returnGraph=TRUE!"))
    
    stateGraph <- stateGraph$graph
  }
    
  geneCols <- setdiff(colnames(stateGraph),c("attractorAssignment","transitionsToAttractor"))
  numGenes <- (length(geneCols)) / 2
  
  from <- apply(stateGraph[,1:numGenes,drop=FALSE],1,paste,collapse="")
  to <- apply(stateGraph[,((numGenes+1):(2*numGenes)),drop=FALSE],1,paste,collapse="")
  vertices <- unique(c(from,to))
  edges <- data.frame(from,to)
  
  res <- graph.data.frame(edges,vertices=as.data.frame(vertices),directed=TRUE)
  res <- set.vertex.attribute(res,"name",value=vertices)    
  
  if ("attractorAssignment" %in% colnames(stateGraph))
    attractorAssignment <- stateGraph$attractorAssignment
  else
  {
    attractorAssignment <- c()
    colorBasins <- FALSE
    drawLegend <- FALSE
  }
  
  if ("stepsToAttractor" %in% colnames(stateGraph))
    attractorIndices <- to[stateGraph$stepsToAttractor == 0]
  else
  {
    attractorIndices <- c()
    highlightAttractors <- FALSE
  }
     
  # determine nodes and edges that belong to attractors
  
  
  # set default edge width and line type
  res <- set.edge.attribute(res,"width",value=0.9)
  res <- set.edge.attribute(res,"lty",value=basin.lty)
  
  if (highlightAttractors)
  {
    attractorEdgeIndices <- which(apply(edges,1,
             function(edge)((edge[1] %in% attractorIndices) & (edge[2] %in% attractorIndices)))) - bias
    # set different edge width and line type for attractor edges
    res <- set.edge.attribute(res,"width",index=attractorEdgeIndices,value=2)
    res <- set.edge.attribute(res,"lty",index=attractorEdgeIndices,value=attractor.lty)
  }

  if (missing(colorSet))
  {
    # define default colors
    colorSet <- c("blue","green","red","darkgoldenrod","gold","brown","cyan",
        "purple","orange","seagreen","tomato","darkgray","chocolate",
        "maroon","darkgreen","gray12","blue4","cadetblue","darkgoldenrod4",
        "burlywood2")
  }
  
  # check for certain graphical parameters in ... 
  # that have different default values in this plot
  if (is.null(args$vertex.size))
    args$vertex.size <- 2
    
  if (is.null(args$edge.arrow.mode))
    args$edge.arrow.mode <- 2
  
   if (is.null(args$edge.arrow.size))
    args$edge.arrow.size <- 0.3
  if (is.null(args$vertex.label.cex))
    args$vertex.label.cex <- 0.5
    
  if (is.null(args$vertex.label.dist))
    args$vertex.label.dist <- 1

  attractors <- unique(attractorAssignment)
  attractors <- attractors[!is.na(attractors)]
  
  if (colorBasins)
  {  
    res <- set.edge.attribute(res,"color",value="darkgrey")
    for (attractor in attractors)
    {
      # determine nodes and edges belonging to the basin of <attractor>
      basinIndices <- which(attractorAssignment == attractor)
      
      if (!symbolic)
      {
        # change vertex color
        res <- set.vertex.attribute(res,"color",basinIndices - bias,
                  value=colorSet[(attractor-1) %% length(colorSet) + 1])
        if (drawLabels)
          res <- set.vertex.attribute(res,"label.color",basinIndices - bias,
              value=colorSet[(attractor-1) %% length(colorSet) + 1])
      }
      #basinEdgeIndices <- which(apply(edges,1,
      #           function(edge)((edge[1] %in% basinIndices) 
      #                   & (edge[2] %in% basinIndices)))) - bias
      basinEdgeIndices <- basinIndices - bias
      
      # change edge color
      res <- set.edge.attribute(res,"color",index=basinEdgeIndices,
              value=colorSet[(attractor-1) %% length(colorSet) + 1])
    }
  }

  if(plotIt)
  {   
    if (drawLabels)
      labels <- vertices
    else
      labels <- NA
    if (piecewise)
      layout <- piecewise.layout(res, layout)
    
    if (symbolic)
      autocurve.edges(res)
    
    do.call("plot",c(list(res),args,"vertex.label"=list(labels), "layout"=list(layout)))
    #plot(res,vertex.size=args$vertex.size,layout=layout,
    #     edge.arrow.mode=args$edge.arrow.mode,
    #     vertex.label=labels,vertex.label.cex=args$vertex.label.cex,
    #     vertex.label.dist=args$vertex.label.dist,
    #     ...)
    if (colorBasins & drawLegend)
      legend(x="bottomleft",pch=15,ncol=1,

             col=colorSet[attractors-1 %% length(colorSet) + 1],
             legend = paste("Attractor",seq_along(attractors)),
             cex=0.5)
  }
  return(invisible(res))
}
