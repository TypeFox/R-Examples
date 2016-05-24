#' Visualize: Generate a graph which vividly displays the gene X, Y and W.
#' 
#' \code{visualize()} generates a graph. It is used to intuitively and vividly display the layout of gene X, Y and W.
#' 
#' @param graph The igraph object of gene network.
#' @param kernel.result The result of graph.kd which finds genes W of a gene X.
#' @param x The gene the plot is generated for.
#' @param k The degree of the neighborhood of X.
#' @param cutoff A threshold to filter gene W.
#' @param path The path where the result graph is saved to. The default path is the original path of input graph.
#' @return a graph of gene X, Y and W
#' @examples \dontrun{
#' kernel <- graph.kd(relate.matrix,g,smoothing.normalize = "one")
#' visualize(g,kernel,x,k=2,cutoff=1,path= NULL)}
#' @export 
#' @importFrom intergraph asNetwork
#' @importFrom GGally ggnet
#' @importFrom ggplot2 ggsave
#' @importFrom stats setNames
visualize <- function(graph, kernel.result, x, k = 2, cutoff = 1, path = NULL) {
  
  X <- as.character(x)
  Y <- V(graph)$name[unlist(igraph::neighborhood(graph, k, nodes = X))]
  z <- kernel.result[X, ]
  W <- names(z[z > cutoff])
  subg <- induced.subgraph(graph, unique(c(X, Y, W)))
  
  type <- vector(mode = "character", length = length(V(subg)))
  type <- setNames(type, V(subg)$name)
  
  for (v in W) {
    type[v] <- "W"
  }
  
  for (v in Y) {
    type[v] <- "Y"
  }
  type[X] <- "X"
  
  network <- asNetwork(subg)
  size <- length(V(subg))
  scale <- (17/961) * size + 1999/961
  
  output <- ggnet(network, node.group = type, segment.size = 1, label.nodes = T, col = "black")
  ggsave(output, filename = paste(as.character(x), ".jpg", sep = ""), path = path, width = 4, height = 3, scale = scale, 
         limitsize = FALSE)
  return(output)
}

#' visualize ego gene X, its k step neighbours, and the W gene communities: Generate a graph 
#' with different community in different colors.
#' \code{visualize.community()}is used to create a graph to display the layout of genes X, X's k-step neighborhood, W and their corresponding community.
#' @param graph The igraph object of gene network.
#' @param kernel.result The result of graph.kd which finds genes W of a gene X. 
#' @param x The gene the plot is generated for.
#' @param k The degree of the neighborhood of X.
#' @param cutoff A threshold to filter gene W.
#' @param community.min The minimum size of the community of W.
#' @param path The path where the result graph is saved to. The default path is the original path of input graph.
#' @return a graph displays genes X, X's k-step neighborhood, and W gene communities in different colors.
#' @examples \dontrun{
#' kernel <- graph.kd(relate.matrix,g,smoothing.normalize = "one")
#' visualize(g,kernel,x,k=2,cutoff=1,community.min=5,path=NULL)}
#' @export
#'
visualize.community <- function(graph, kernel.result, x, k = 2, cutoff = 1, community.min = 5, path = NULL) {
  X <- as.character(x)
  Y <- V(graph)$name[unlist(igraph::neighborhood(graph, k, nodes = X))]
  z <- kernel.result[X, ]
  
  W <- names(z[z > cutoff])
  
  wc <- getCommunity(z, graph, cutoff, community.min)
  member <- membership(wc)
  community_index <- names(sizes(wc)[sizes(wc) > community.min])
  if (length(community_index) > 7) {
    print("comunity too large.")
    return()
  }
  
  subg <- induced.subgraph(graph, unique(c(X, Y, W)))
  
  type <- rep("other", length(V(subg)))
  type <- setNames(type, V(subg)$name)
  
  time <- 1
  for (ci in community_index) {
    w <- names(member[member == ci])
    type[w] = paste("w", time, sep = "")
    time <- time + 1
  }
  
  type[Y] <- "Y"
  
  type[X] <- "X"
  
  network <- asNetwork(subg)
  size <- length(V(subg))
  scale <- (17/961) * size + 1999/961
  
  output <- ggnet(network, node.group = type, segment.size = 1, label.nodes = T, col = "black")
  
  ggsave(output, filename = paste(as.character(x), ".jpg", sep = ""), path = path, width = 4, height = 3, scale = scale, 
         limitsize = FALSE)
  return(output)
}