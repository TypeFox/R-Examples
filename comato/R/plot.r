min.distance <- function(centers, com, layouts)
{
  centers = matrix(centers, ncol=2)
  layout = cbind(rep(0, length(com$membership)), rep(0, length(com$membership)))
  
  for (i in 1:max(com$membership))
  {
    nodes <- which(com$membership == i)     
    layout[nodes,1] =  layouts[[i]][,1] + centers[i,1]
    layout[nodes,2] =  layouts[[i]][,2] + centers[i,2]                     
  }  
  return(min(dist(layout, method="manhattan")))
}

flip.blocks <- function(centers, com, layouts)
{
  centers = matrix(centers, ncol=2)
  cp = sample(max(com$membership)) 
  return(centers[cp,]) 
}
  
calculate.degrees <- function(nodeCount) {
  return(0:(nodeCount-1)*2*pi/nodeCount)
}

calculate.layout <- function(map)
{
  layout = cbind(rep(0, length(V(map))), rep(0, length(V(map))))

  com <- walktrap.community(map)  
  
  degs = calculate.degrees(max(com$membership)-1)
  centers = log2(2*max(summary(com$membership)))*cbind(cos(degs), sin(degs))
  centers = rbind(c(0, 0), centers)
  layouts = list()
  
  for (i in 1:max(com$membership)) 
  {  
    nodes <- which(com$membership == i)		
    res <- layout.fruchterman.reingold(induced.subgraph(map, nodes))
    res[,1] = log2(length(nodes))*((res[,1] - min(res[,1])) / (max(res[,1]) - min(res[,1])) - 1)
    res[,2] = log2(length(nodes))*((res[,2] - min(res[,2])) / (max(res[,2]) - min(res[,2])) - 1)
    
    layouts = c(layouts, list(res))   
  }

  if (max(com$membership) > 4)
  {
    centers = optim(centers, min.distance, flip.blocks, com=com, layouts=layouts, method="SANN", control = list(fnscale=-1))$par 
    centers = matrix(centers, ncol=2)    
  }
  
  for (i in 1:max(com$membership))
  {
    nodes <- which(com$membership == i)     
    layout[nodes,1] =  layouts[[i]][,1] + centers[i,1]
    layout[nodes,2] =  layouts[[i]][,2] + centers[i,2]                     
  }  
  
  return(layout)
}


#' Plotting a conceptmap
#' 
#' \code{plot} plots a concept map. Includes finding a good layout based on communities and a circular layout.
#' Is especially suited for plotting larger concept maps, in particular amalgamations.
#' @param x A conceptmap object.
#' @param edge.labels If TRUE, the labels of edges will be plotted as well.
#' @param max.label.len The maximal length of labels (in characters) that are plotted completely. Longer labels will be shortend by "...".
#' @param scale Overall scaling factor that is applied to the plot.
#' @param layout If not NULL, must either be one of "fruchterman.reingold", "kamada.kawai", "spring" or "reingold.tilford" or a numeric matrix.
#' If it is a string, the corresponding layouting algorithm of the igraph package will be called. If it is a numeric matrix, it must contain a row
#' for each concept and two columns that determine the x and y coordinates of this concept. Then this layout will be used directly.
#' @param ... -
#' @return -
#' @examples
#' #Create concept map from a random graph
#' require("igraph")
#' g1 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' E(g1)$name <- rep("", length(E(g1)))
#' plot(conceptmap(g1), edge.labels=FALSE, layout="kamada.kawai")
#' @method plot conceptmap
#'@export 
plot.conceptmap <- function(x, edge.labels=T, max.label.len=25, scale=1, layout=NULL, ...)
{
  t = x$map
  temp <- substr(V(t)$name, 1, max.label.len)
  temp[which(nchar(V(t)$name) > max.label.len)] = paste( temp[which(nchar(V(t)$name) > max.label.len)], "...", sep="")
  V(t)$label <- temp
  E(t)$label <- E(t)$name  
  if (!edge.labels)
    t = remove.edge.attribute(t, "label")
  V(t)$size <- 8*scale + 2.75*scale*nchar(temp)
  if (is.null(layout))
    l = calculate.layout(t)
  else
    if (is.character(layout))
    {
      if (layout == "fruchterman.reingold")
        l = layout.fruchterman.reingold(t)  
      if (layout == "kamada.kawai")
        l = layout.kamada.kawai(t)        
      if (layout == "spring")
        l = layout.spring(t)        
      if (layout == "reingold.tilford")
        l = layout.reingold.tilford(t)              
    }
  else
    l = layout
  old.p <- par("mai")
  par(mai=c(0,0,0,0))
  plot(t, layout=l, 
       vertex.color= "#dad7cb", vertex.shape = "crectangle", vertex.size2=8*scale, vertex.label.color="black", vertex.frame.color="#dad7cb", 
       edge.color="black", edge.arrow.size=0.25*scale, vertex.label.cex=0.55*scale)
  par(mai=old.p)
}


#' Plotting a series of concept maps
#' 
#' \code{plot} plots a set of concept maps. The layout is determined based on the union of all concept maps, then each
#' map is individually plotted using this fixed layout. Is escpecially useful for visualizing horizontal landscapes.
#' @param x A conceptmaps object.
#' @param edge.labels If TRUE, the labels of edges will be plotted as well.
#' @param max.label.len The maximal length of labels (in characters) that are plotted completely. Longer labels will be shortend by "...".
#' @param scale Overall scaling factor that is applied to the plot.
#' @param layout If not NULL, must be one of "fruchterman.reingold", "kamada.kawai", "spring" or "reingold.tilford".
#' The corresponding layouting algorithm of the igraph package will be called. If it is NULL, the layouting based on communities and
#' a circular layout will be used.
#' @param ... -
#' @return -
#' @examples
#' #Create concept maps from three random graphs
#' require("igraph")
#' g1 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' g2 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' g3 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' E(g1)$name <- rep("", length(E(g1)))
#' E(g2)$name <- rep("", length(E(g2)))
#' E(g3)$name <- rep("", length(E(g3)))
#' #Create conceptmaps object from three conceptmap objects
#' simple_cms = conceptmaps(list(conceptmap(g1), conceptmap(g2), conceptmap(g3)))
#' 
#' plot(simple_cms, layout="spring")
#' @method plot conceptmaps
#'@export 
plot.conceptmaps <- function(x, edge.labels=T, max.label.len=25, scale=1, layout=NULL, ...)
{
  t = sort(get.unified.concepts(x))
  maps = unify.concepts(x)
  graphs = list()
  for (m in maps$maps)
  {
    adj = get.adjacency(m$map)
    i = sort.int(colnames(adj), index.return=T)$ix
    adj = adj[i, i]
    graphs = c(graphs, list(graph.adjacency(adj)))    
  } 
  union = graph.union(graphs)
  V(union)$name <- t
  V(union)$label <- t
  if (is.character(layout))
  {
    if (layout == "fruchterman.reingold")
      l = layout.fruchterman.reingold(union)  
    if (layout == "kamada.kawai")
      l = layout.kamada.kawai(union)        
    if (layout == "spring")
      l = layout.spring(union)        
    if (layout == "reingold.tilford")
      l = layout.reingold.tilford(union)              
  }
  else
    l = calculate.layout(union)
  for (m in maps$maps)
  {
    concepts = match(V(m$map)$name, t)
    plot(m, layout=l[concepts,], edge.labels=edge.labels, scale=scale, max.label.len=25)    
  }
}


