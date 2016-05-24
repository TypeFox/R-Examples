#' Analyzing graph measures of a concept map
#' 
#' \code{analyze.graph.measures} analyzes several basic graph measures of a given graph in form of a conceptmap object. All measures
#' are derived by the appropriate functions of igraph.
#' @param x A conceptmap object.
#' @return A list with named components that contain the betweenness measure, the edge.connectivity, the diameter, the degree distribution and the 
#' communities using the Fastgreedy algorithm.
#' @examples 
#' require("igraph")
#' g1 = set.vertex.attribute(erdos.renyi.game(15, 0.7, type="gnp"), "name", value=1:15)
#' analyze.graph.measures(conceptmap(g1))
#' @method analyze.graph.measures conceptmap
#'@export 
analyze.graph.measures.conceptmap <- function(x) 
{
  return(analyze.graph.measures(x$map))
}


#' Analyzing graph measures of an igraph object
#' 
#' \code{analyze.graph.measures.igraph} is a convenience method that can be called directly on the result of \code{\link{landscape}}.
#' It works just like \code{\link{analyze.graph.measures.conceptmap}}
#' @param x An igraph object.
#' @return A list with named components that contain the betweenness measure, the edge.connectivity, the diameter, the degree distribution and the 
#' communities using the Fastgreedy algorithm.
#' @examples 
#' require("igraph")
#' g1 = set.vertex.attribute(erdos.renyi.game(15, 0.7, type="gnp"), "name", value=1:15)
#' analyze.graph.measures(g1)
#' @method analyze.graph.measures igraph
#'@export 
analyze.graph.measures.igraph <- function(x) 
{
  res <- list(betweenness(x), edge.connectivity(x), diameter(x), degree(x), fastgreedy.community(as.undirected(x, mode="collapse")))
  names(res) <- c("betweenness", "connectivitiy", "diameter", "degrees", "communities")
  return(res)
}


#' Analyzing small subgraph patterns.
#' 
#' \code{analyze.subgraphs} analyzes the frequency of subgraph patterns given a list of concepts and a set of maps in a conceptmaps object.
#' @param maps A conceptmaps object.
#' @param concept.list A list of concepts (as strings) that define the subgraphs to be analyzed. must be between 2 and 4 concepts in length.
#' @return A list with two elements. The first is a vector that contains the absolute number of occurrences for each subgraph pattern.
#' The second element is a list of igraph objects of the pattern themselves. 
#' @examples 
#' #Create concept maps from three random graphs
#' require("igraph")
#' g1 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' g2 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' g3 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' 
#' #Create conceptmaps object from three conceptmap objects
#' simple_cms = conceptmaps(list(conceptmap(g1), conceptmap(g2), conceptmap(g3)))
#' 
#' analyze.subgraphs(simple_cms, c("1", "2", "3"))
#'@export 
analyze.subgraphs <- function(maps, concept.list) 
{
  stopifnot(length(concept.list) > 1 || length(concept.list) < 5)
  temp = modify.concepts(maps, concept.list)
  n = length(concept.list)
  graphs = c()
  t = matrix(0, nrow=n, ncol=n, dimnames=list(concept.list, concept.list))
  pos = length(which(upper.tri(t)))
  for (i in 1:2^((n^2-n)/2))
  {
    t[upper.tri(t)] <-  as.numeric(intToBits(i))[1:pos]
    graphs= c(graphs, list(graph.adjacency(t, mode="upper")))  
  }
  
  res = rep(0, 2^((n^2-n)/2))
  for (m in temp$maps) {
    a = get.adjacency(as.undirected(simplify(m$map)))
    ix_a = sort.int(colnames(a), index.return=T)$ix
    for (g in 1:length(graphs)) {
      b = get.adjacency(graphs[[g]])
      ix_b = sort.int(colnames(b), index.return=T)$ix
      if (all(a[ix_a, ix_a] == b[ix_b, ix_b]))
      {
        res[g] = res[g] + 1
        break
      }
    }
  }
  return(list(res, graphs))
}


#' Analyzing graph similarity.
#' 
#' \code{analyze.similarity} calculates a measure of graph similarity between two concept maps.
#' @param map1 A conceptmap object.
#' @param map2 A conceptmap object.
#' @return A value between 0 and 1 that indicated the structural similariy of the underlying graphs. The graphs need not 
#' share the same set of nodes.
#' @examples 
#' require("igraph")
#' g1 = set.vertex.attribute(erdos.renyi.game(15, 0.7, type="gnp"), "name", value=1:15)
#' g2 = set.vertex.attribute(erdos.renyi.game(15, 0.7, type="gnp"), "name", value=1:15)
#' analyze.similarity(conceptmap(g1), conceptmap(g2))
#' @seealso The structural similarity that is calculated is described in: 
#' Goldsmith, Timothy E.; Davenport, Daniel M. (1990): Assessing Structural Similarity of Graphs. In: Roger W. Schvaneveldt (Hg.): Pathfinder associative networks. Studies in knowledge organizations. Norwood, N.J: Ablex Pub. Corp., S. 74-87.
#'@export 
analyze.similarity <- function(map1, map2)
{
  m1 = get.adjacency(as.undirected(map1$map, mode="collapse"), type="both")
  m2 = get.adjacency(as.undirected(map2$map, mode="collapse"), type="both")
  e = intersect(colnames(m1)[which(colSums(m1) == 0)], colnames(m2)[which(colSums(m2) == 0)])
  sum = length(e)
  for (v in setdiff(colnames(m1), e))
  {
    n1 = names(which(m1[v,] > 0))
    n2 = names(which(m2[v,] > 0))
    sum = sum + length(intersect(n1, n2))/length(union(n1, n2))
  }
  return(sum/dim(m1)[1]) 
}
