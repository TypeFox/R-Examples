# Plots a graph of the transitions in a probabilistic Boolean network.
# <markovSimulation> is the result of a Markov simulation with returnTable=TRUE.
# If stateSubset is specified, only transitions between states in the set are considered.
# If <drawProbabilities> is true, the edges are annotated with probabilities.
# If <drawStateLabels> is true, the states are annotated with their gene values.
# <layout> specifies the igraph layout to be used.
# If <plotIt> is false, only the graph is returned, and nothing is plotted.
# ... specifies further parameters to igraph.
plotPBNTransitions <- function(markovSimulation,stateSubset,
                               drawProbabilities=TRUE,drawStateLabels=TRUE,
                               layout=layout.fruchterman.reingold,
                               plotIt=TRUE,...)
{
  stopifnot(inherits(markovSimulation,"MarkovSimulation"))

  if (is.null(markovSimulation$table))
    stop(paste("The supplied simulation result does not contain transition information.",
               "Please re-run markovSimulation() with returnTable=TRUE!"))
  
  # assemble edges from table
  edgeMatrix <- data.frame(apply(markovSimulation$table$initialStates,2,
                           function(x)paste(dec2bin(x,length(markovSimulation$genes)),collapse="")),
                           apply(markovSimulation$table$nextStates,2,
                           function(x)paste(dec2bin(x,length(markovSimulation$genes)),collapse="")))
                           
  if (!missing(stateSubset))
  {    
    # determine edges to be excluded based on the subset
    stateSubset <- sapply(stateSubset,function(x)paste(x,collapse=""))
    keepIndices <- apply(edgeMatrix,1,function(row)
                        {
                          (length(intersect(row,stateSubset)) == length(unique(row)))
                        })
    
    # drop edges                   
    edgeMatrix <- edgeMatrix[keepIndices,]
    probabilities <- markovSimulation$table$probabilities[keepIndices]
  }
  else
    probabilities <- markovSimulation$table$probabilities
  
  # determine set of vertices
  vertices <- as.data.frame(as.character(unique(c(as.character(edgeMatrix[,1]),
                            as.character(edgeMatrix[,2])))))

  # build graph                             
  graph <- graph.data.frame(edgeMatrix,vertices=vertices,directed=TRUE)
  
  if (drawProbabilities)
    graph <- set.edge.attribute(graph,"label",value=paste("    ",probabilities))
  if (drawStateLabels)
    label <- as.character(vertices[,1])
  else
    label <- NA
    
    
  if (plotIt)
  {
    # set default values for further graphical parameters
    args <- list(...)
  
    if (is.null(args$vertex.size))
      args$vertex.size <- 2

    if (is.null(args$edge.arrow.mode))
      args$edge.arrow.mode <- 0

    if (is.null(args$vertex.label.cex))
      args$vertex.label.cex <- 0.75
      
    if (is.null(args$edge.label.cex))
      args$edge.label.cex <- 0.75      

    if (is.null(args$vertex.label.dist))
      args$vertex.label.dist <- 1 

    if (is.null(args$vertex.color))
      args$vertex.color <- "grey"

    if (is.null(args$edge.label.color))
      args$edge.label.color <- "green"

    if (is.null(args$edge.arrow.size))
      args$edge.arrow.size <- 0.5 

    # plot it       
    plot(graph,vertex.label=label,layout=layout,vertex.label.cex=args$vertex.label.cex,
         vertex.size=args$vertex.size,vertex.color=args$vertex.color,
         vertex.label.dist = args$vertex.label.dist,
         edge.arrow.size=args$edge.arrow.size,
         edge.label.color=args$edge.label.color,
         edge.label.cex=args$edge.label.cex,...)
  }
  
  # return the igraph object
  return(invisible(graph))
}

