igraph.to.gexf <- function(igraph.obj, position=NULL) {
  
  g <- igraph.obj
  rm(igraph.obj)
  
  # Retrive elements from igraph object
  gdata <- get.data.frame(g, what="both")
  tmpedges <- gdata$edges
  tmpnodes <- gdata$vertices
  
  # If nodes have no name
  if (!("name" %in% colnames(tmpnodes))) 
    tmpnodes <- data.frame(tmpnodes, name=1:nrow(tmpnodes))
  
  # Nodes and edges list
  nodes <- edge.list(tmpedges[,c(1,2)])
  edges <- nodes$edges
  nodes <- nodes$nodes
  
  # Building nodes
  if (length(tmpnodes)) {
    nodes <-merge(tmpnodes,nodes, by.x="name", by.y="label")
#    nodes <- nodes[,!(colnames(nodes) %in% c("label"))]
  }
  
  # Nodes Attributes
  x <- list.vertex.attributes(g)
  x <- x[!(x %in% c("label","color","size"))]
  if (!length(x)) nAtt <- NULL
  else nAtt <- subset(nodes, select=x)
  
  # Edges Attributes
  x <- list.edge.attributes(g)
  x <- x[!(x %in% c("weight","color","edgesLabel","width"))]
  if (!length(x)) eAtt <- NULL
  else eAtt <- subset(tmpedges, select=x)
  
  # Edges Weights
  if (length(E(g)$weight)) eW <- E(g)$weight
  else eW <- NULL
  
  # Nodes Viz att
  if (length(tmpnodes$color)) {
    nVizAtt <- list(color=t(col2rgb(tmpnodes$color, alpha=T)))
  }
  else nVizAtt <- NULL
  
  if (length(tmpnodes$size)){
    nVizAtt$size <- tmpnodes$size
  }
  
  nVizAtt$position <- position
  
  # Edges Viz att
  if (length(tmpedges$color)) {
    eVizAtt <- list(color=t(col2rgb(tmpedges$color, alpha=T)))
  }
  else eVizAtt <- NULL
  
  # Edge type
  if (is.directed(g)) defaultedgetype <- "directed"
  else defaultedgetype <- "undirected"
  
  # Building graph
  return(
    write.gexf(
      nodes = nodes[,c("id","name")], 
      edges = edges, 
      edgesAtt = eAtt,
      nodesAtt = nAtt,
      defaultedgetype = defaultedgetype,
      nodesVizAtt = nVizAtt,
      edgesVizAtt = eVizAtt,
      edgesWeight = eW
      )
    )
}

gexf.to.igraph <- function(gexf.obj) {

  # Checks the class
  if (!inherits(gexf.obj,"gexf")) stop("-graph- is not of -gexf- class.") 
 
  g <- gexf.obj
  rm(gexf.obj)
  #colnames(g$nodes)[colnames(g$nodes)=="label"] <- "name"
  
  # Starting igraph object
  #colnames(g$nodes)[colnames(g$nodes) == "name"]
  
  g2 <- graph.data.frame(
    g$edges[,c("source","target")],
    directed=(g$mode[[1]] != "undirected"),
    vertices=g$nodes[,c("id"),drop=FALSE]
    )
  
  # Labels
  V(g2)$name <- g$nodes$label
  
  # Nodes Viz atts
  if (length(x <- g$nodesVizAtt$color)) {
    V(g2)$color <- rgb(x$r/255,x$g/255,x$b/255,x$a/255)
  }

  # Nodes atts
  if (length(x <- g$nodes[, !(colnames(g$nodes) %in% c("id", "label"))])) {
    for(i in names(x)) g2 <- set.vertex.attribute(g2, i, value=x[,c(i)])
  }

  # Edges atts
  if (length(x <- g$edges[, !(colnames(g$edges) %in% c("id", "source", "target", "weight") )])) {
    for(i in names(x)) g2 <- set.edge.attribute(g2, i, value=x[,c(i)])
  } 
  
  # Edges Viz atts
  if (length(x <- g$edgesVizAtt$color)) {
    E(g2)$color <- rgb(x$r/255,x$g/255,x$b/255,x$a/255)
  }
  
  # Edges weights 
  if (length(x <- g$edges$weight)) {
   E(g2)$weight <- x
  }
  
  return(g2)
}

