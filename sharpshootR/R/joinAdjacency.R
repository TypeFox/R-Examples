joinAdjacency <- function(x, vars=c('l_musym', 'r_musym')) {
  
  # extract "left" and "right" map unit symbols, removing missing values
  d <- slot(x, 'data')[, vars]
  edge.list <- as.matrix(na.omit(d))
  
  # init igraph object: note that there will be many duplicate edges
  g <- graph.edgelist(edge.list, directed=FALSE)
  
  # keep track of duplicate edges as weight, then remove
  E(g)$weight <- log(count.multiple(g))
  g <- simplify(g)
  
  # save as weighted adjacancy matrix for plotting with sharpshootR functions
  a <- get.adjacency(g, attr='weight')
  
  return(a)
}