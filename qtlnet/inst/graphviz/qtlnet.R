## This was written to handle graphviz for qtlnet routines but was later abandoned.
## Can this be rewritten to use graph.qtlnet output? (to reduce code duplication)

Rgraphviz.qtlnet <- function(graph, layout.method = "dot", ...)
{
  require(Rgraphviz)
  
  mygR <- graph.and.attributes(graph, ...)

  ## Plot the graph object.
  plot(mygR$graph, layout.method, edgeAttrs = mygR$edges, nodeAttrs = mygR$nodes,
                 attrs = mygR$all.nodes)
}
############################################################
graph.and.attributes <- function(graph,
                                 node.shape = "ellipse",
                                 lwd = 3, fontsize = 18,
                                 width = 3.0, height = 1.0, ...)
{

  gR <- create.directed.graph.object(graph)

  ## Edge attributes.
  graph.edges <- get.graph.edges(graph)
  n.edges <- nrow(graph.edges)
  edges.names <- paste(graph.edges$cause, graph.edges$effect, sep = "~")
  tmp <- pmin(20, 1 + floor(20 * 3 * (1 - graph.edges$width) / 2))
  aux <- grey(tmp / 20)
##  aux[auxDG1$width < 1/3] <- NA
  eAttrs <- list(color = aux,
                 lwd = as.character(rep(lwd, n.edges)))
  for(i in names(eAttrs))
    names(eAttrs[[i]]) <- edges.names

  ## Node attributes.
  vertices <- get.graph.vertices(graph)
  n.vert <- nrow(vertices)
  nAttrs <- list(fillcolor = as.character(vertices$color),
                 fontsize = as.character(rep(fontsize, n.vert)),
                 width = as.character(rep(width, n.vert)),
                 height = as.character(rep(height, n.vert)))
  for(i in names(nAttrs))
    names(nAttrs[[i]]) <- vertices$name

  ## Shape attribute.
  attrs <- list(node =list(shape = node.shape))

  list(graph = gR, edges = eAttrs, nodes = nAttrs,
       all.nodes = attrs)
}
################################################################
create.directed.graph.object <- function(graph)
{
  
  mynodes <- as.character(get.graph.vertices(graph)$name)
  output <- get.graph.edges(graph)
  n.edges <- nrow(output)
  le <- length(mynodes)
  edL <- vector("list",length=le)
  names(edL) <- mynodes

  ## Can probably improve on this, but it is low cost.
  for(i in 1:le){
    auxNode <- mynodes[i]
    auxEdges <- c()
    for(j in 1:n.edges) {
      if(output[j,1] == auxNode){
        auxEdges <- c(auxEdges,which(mynodes==output[j,2]))
      }
    }
    edL[[i]] <- list(edges = auxEdges)
  }
  new("graphNEL", nodes = mynodes, edgeL = edL, edgemode = "directed")
}
