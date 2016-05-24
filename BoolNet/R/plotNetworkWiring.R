# Plot a wiring graph of the network <network> with the supplied
# graphical parameters.
# Requires igraph.
# Returns the igraph structure representing the wiring graph. 
plotNetworkWiring <- function(network,layout=layout.fruchterman.reingold,plotIt=TRUE,...)
{
  stopifnot(inherits(network,"ProbabilisticBooleanNetwork") | inherits(network,"BooleanNetworkCollection")
            | inherits(network,"BooleanNetwork") | inherits(network,"SymbolicBooleanNetwork"))
  
  if (installed.packages()["igraph","Version"] < package_version("0.6"))
    bias <- 1
  else
    bias <- 0 
    
  edgeList <- c()
  
  # construct list of edges from interactions

  if (inherits(network,"BooleanNetwork"))
  # deterministic network
  {
    for (i in seq_along(network$genes))
    {
      if (network$interactions[[i]]$input[1] != 0)
      # no edges for constant genes
      {
        edgeList <- rbind(edgeList,
                  cbind(network$interactions[[i]]$input,
                  rep(i,length(network$interactions[[i]]$input))))
      }
    }
  }
  else
  if (inherits(network,"SymbolicBooleanNetwork"))
  # symbolic network
  {
    inputs <- lapply(network$interactions, getInputs, index=TRUE)
    for (i in seq_along(network$genes))
    {
      edgeList <- rbind(edgeList,
                  cbind(inputs[[i]],
                  rep(i,length(inputs[[i]]))))
    }
  }  
  else
  # probabilistic network
  {
    for (i in seq_along(network$genes))
    {
      for (j in seq_along(network$interactions[[i]]))
      {
        if (network$interactions[[i]][[j]]$input[1] != 0)
        # no edges for constant genes
        {
          edgeList <- rbind(edgeList,
                            cbind(network$interactions[[i]][[j]]$input,
                            rep(i,length(network$interactions[[i]][[j]]$input))))
        }
      }
    }
  }

  # build graph from edge list
  res <- graph.data.frame(edgeList-bias,directed=TRUE,vertices=as.data.frame((seq_along(network$genes)) - bias))
  res <- set.vertex.attribute(res,"name",value=network$genes)
  
  args <- list(...)
  
  # check for certain graphical parameters in ... 
  # that have different default values in this plot
  if (is.null(args$vertex.color))
    args$vertex.color <- "grey"
    
  if (is.null(args$edge.arrow.size))
    args$edge.arrow.size <- 0.5
    
  if (is.null(args$vertex.label.cex))
    args$vertex.label.cex <- 0.7

  if (is.null(args$vertex.size))
    args$vertex.size <- 18    

  if (plotIt)
  {
    plot(res,vertex.label=network$genes,vertex.label.cex=args$vertex.label.cex,
         vertex.size=args$vertex.size,vertex.color=args$vertex.color,
         edge.arrow.size=args$edge.arrow.size,
         layout=layout,...)
  }
  return(invisible(res))
}
