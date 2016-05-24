#' Forming the concept vector of a conceptmap object
#' 
#' \code{concept.vector} transforms a concept map into a numeric vector that contains for each occuring concept the number of 
#' adjacent edges.
#' @param x A conceptmap object.
#' @return A numeric vector. The columns are named after the concepts and sorted alphabetically.
#' @examples
#' #Create concept map from a random graph
#' require("igraph")
#' g1 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' cm = conceptmap(g1)
#' 
#' concept.vector(cm)
#'@export 
concept.vector <- function(x) 
{
  res = as.numeric((rowSums(get.adjacency(x$map)) + colSums(get.adjacency(x$map))) > 0)
  names(res) <- x$concepts
  ix = sort.int(names(res), index.return=T)$ix
  return(res[ix])
}


#' Forming the edge vector of a conceptmap object
#' 
#' \code{edge.vector} transforms a concept map into a numeric vector that contains for each occuring pair of concepts whether or not
#' this pair is connected by a proposition in the concept map.
#' @param x A conceptmap object.
#' @return A numeric vector. The columns are named after the concept-pairs which are sorted alphabetically.
#' @examples
#' #Create concept map from a random graph
#' require("igraph")
#' g1 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' cm = conceptmap(g1)
#' 
#' edge.vector(cm)
#'@export 
edge.vector <- function(x) 
{
  adj = get.adjacency(x$map)
  ix = sort.int(x$concepts, index.return=T)$ix
  res = as.vector(adj[ix, ix])
  temp = cbind(rep(sort(x$concepts), length(x$concepts)), rep(sort(x$concepts), times=1, each=length(x$concepts)))
  names(res) <- apply(temp, MARGIN=1, FUN=paste, collapse="-")
  return(res)
}


#' Aggregating the maps of a conceptmaps object into a concept landscape
#' 
#' \code{landscape} transforms a set of concept maps into a concept landscape using one of several possible methods.
#' Depending on the value of \code{result} and accumulation or amalgamation is performed on the concept map data. The amalgamation
#' forms a weighted graph based on the unified set of concepts. An accumulation transforms each concept map into a vector and
#' returns a matrix of these vectors. Using \code{FUN} the process of transformation can be influenced in both cases.
#' @param maps A conceptmaps object.
#' @param result Either "graph" or "matrix" defines the return type and the method of aggregation. An amalgamation results in a weighted graph and
#' an accumulation results a matrix.
#' @param mode If \code{result} is "graph", it can be either "directed" or "undirected" deciding how the graph should be formed. 
#' First, a weight matrix is formed from the set of concept maps. If \code{FUN} is NULL, the weights simply reflect the number of maps in which 
#' a given edge is present. Otherwise, \code{FUN} must be a function that accepts three parameters and returns a numeric value. For the first parameter
#' the current conceptmap object will be passed and for the second and third parameters the start and end concepts of an edge is passed.
#' The return value of the function will then be added to the weight matrix. If \code{result} is "matrix", it can be one of "graph.sim", "concept.vector", "edge.vector" or "custom". "graph.sim" return the graph similarity matrix,
#' "concept.vector" and "edge.vector" return the concept matrix or edge matrix respectively. Finally, "custom" can be used for arbitrary
#' transformations: For each map, the function passed to \code{FUN} and the resulting vector is forming a row of the returned matrix.
#' @param FUN If \code{result} is "matrix" and \code{mode} is "custom", a function with one parameter of type conecptmap must be given that is 
#' called for each of the consituent maps. The function must return a numeric vector of equal length for all maps of a conceptmaps object.
#' @return Depending on \code{result} either an igraph object or a numeric matrix.
#' @examples
#' #Create concept maps from three random graphs
#' require("igraph")
#' g1 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' g2 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' g3 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' 
#' #Create conceptmaps object from three conceptmap objects
#' cms = conceptmaps(list(conceptmap(g1), conceptmap(g2), conceptmap(g3)))
#' 
#' landscape(cms, result="graph", mode="undirected")
#' 
#' landscape(cms, result="matrix", mode="concept.vector")
#'@export 
landscape <- function(maps, result = c("graph", "matrix"), mode, FUN=NULL)
{
  result = match.arg(result)
  if (result == "graph")
  {
    if (is.null(FUN))
      FUN = function(x, e1, e2) {return(1)}
    concepts = get.unified.concepts(maps)
    mat = matrix(0, nrow=length(concepts), ncol=length(concepts), dimnames=list(concepts, concepts))
    for (m in 1:length(maps$maps))
    {
      res = get.edgelist(maps$maps[[m]]$map)
      if (dim(res)[1] > 0)
        for (i in 1:dim(res)[1])
          mat[which(concepts == res[i, 1]), which(concepts == res[i, 2])] = mat[which(concepts == res[i, 1]), which(concepts == res[i, 2])] + FUN(maps$maps[[m]], res[i, 1], res[i, 2])
    }
    if (mode == "undirected")
      return(graph.adjacency(mat, mode="plus", weighted=T))
    else
      return(graph.adjacency(mat, mode="directed", weighted=T))
  }
  else if (result == "matrix")
  {
    temp = unify.concepts(maps)
    if (mode == "graph.sim")
    {
      mat = matrix(0, nrow=temp$count, ncol=temp$count)
      for (i in 1:temp$count)
        for (j in i:temp$count)
          mat[i, j] = analyze.similarity(temp$maps[[i]], temp$maps[[j]])
      mat[lower.tri(mat)] = mat[upper.tri(mat)]
      return(mat)
    }
    mat = c()
    for (m in temp$maps)
    {
      col = c()
      if (mode == "concept.vector")
        col = concept.vector(m)
      else if (mode == "edge.vector")
        col = edge.vector(m)        
      else if (mode == "custom")
        col = FUN(m)        
      mat = rbind(mat, col)
    }
    return(mat)
  }
}
