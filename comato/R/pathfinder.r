#' Creating a Pathfinder network from a matrix
#' 
#' \code{pathfinder} creates the Pathfinder network from a given weight matrix.
#' @param data A non-negative weight matrix of a graph that can be either directed or undirected. 
#' @param q The parameter q used in the Pathfinder algorithm. The resulting graph will be q-triangular.
#' @param r The parameter r used in the Pathfinder algorithm for the r-metric.
#' @param ... -
#' @return A numeric weight matrix that represented the Pathfinder graph of the input graph.
#' @seealso The Pathfinder algorithm is implemented based on the description in:
#' Dearholt, Donald W.; Schvaneveldt, Roger W. (1990): Properties of Pathfinder Netowrks. In: Roger W. Schvaneveldt (Hg.): Pathfinder associative networks. Studies in knowledge organizations. Norwood, N.J: Ablex Pub. Corp., S. 1-30.
#' @examples 
#' #Manually create a weighted graph
#' data = matrix(data = 0, nrow = 6, ncol=6)
#' colnames(data) <- c("Object", "Class", "Method", "Attribute", "Visibility", "Algorithm")
#' rownames(data) <- c("Object", "Class", "Method", "Attribute", "Visibility", "Algorithm")
#' data["Object", "Class"] = 3
#' data["Object", "Method"] = 3
#' data["Object", "Attribute"] = 10
#' data["Object", "Visibility"] = Inf
#' data["Object", "Algorithm"] = 9
#' data["Class", "Method"] = 7
#' data["Class", "Attribute"] = 6
#' data["Class", "Visibility"] = 8
#' data["Class", "Algorithm"] = 10
#' data["Method", "Attribute"] = 4
#' data["Method", "Visibility"] = 9
#' data["Method", "Algorithm"] = 3
#' data["Attribute", "Visibility"] = 5
#' data["Attribute", "Algorithm"] = 10
#' data["Visibility", "Algorithm"] = Inf
#' 
#' data = data + t(data)
#' 
#' #Run the Pathfinder algorithm with several different parameters
#' pathfinder(data, 5, 1)
#' pathfinder(data, 2, 1)
#' pathfinder(data, 5,Inf)
#' pathfinder(data, 2, Inf)
#' @method pathfinder matrix
#'@export
pathfinder.matrix <- function(data, q, r, ...)
{
  n = dim(data)[1]
  W = signif(data)
  D = signif(data)
  for (i in 2:q)
  {
    old_data = W
    for (j in 1:n)
      for (k in 1:n)
      {
        if (r == Inf)
          W[j, k] = signif(min(pmax(data[j,], old_data[,k])))
        else
          W[j, k] = signif(min(data[j,]^r + old_data[,k]^r)^(1/r))
        if (D[j, k] > W[j, k])
          D[j, k] = W[j, k]
      }    
  }
  for (j in 1:n)
    for (k in 1:n)
      if (D[j, k] < signif(data[j, k]))
        W[j, k] = Inf
  return(W)
}


#' Creating a Pathfinder network from a conceptmaps object
#' 
#' \code{pathfinder} creates the Pathfinder network from a given set of conceptmaps. The conecpts of each concept map are
#' unified, then the concept maps are transformed into a weight matrix and \code{\link{pathfinder.matrix}} is called on the data.
#' @param data A conceptmaps object.
#' @param q The parameter q used in the Pathfinder algorithm. The resulting graph will be q-triangular.
#' @param r The parameter r used in the Pathfinder algorithm for the r-metric.
#' @param threshold A numeric value used for pruning the graph before the Pathfinder algorithm. The pruning works in conjunction with
#' the value of \code{prune.edges}.
#' @param prune.edges If TRUE, each entry of the weight matrix that is lower than \code{threshold} will be set to 0 and columns
#' with a resulting sum of 0 are removed. If FALSE, only columns of the weight matrix with a sum of less than \code{threshold} 
#' will be removed.
#' @param return.cm If TRUE, a conceptmap object will be returned. Otherwise, an igrpah object will be returned.
#' @param filename Optional. If specified, the resulting Pathfinder network will be stored in TGF format in the given file.
#' @param ... -
#' @return Depending on \code{return.cm} either an igraph object or a conceptmap object that represents the Pathfinder network.
#' If an igraph object is returned, the graph will be weighted.
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
#' #Create Pathfinder network from data and return a conceptmap object
#' cm = pathfinder(simple_cms, q=1, return.cm=TRUE)
#' @method pathfinder conceptmaps
#'@export
pathfinder.conceptmaps <- function(data, q=Inf, r=1, threshold=0, prune.edges=F, return.cm=F, filename="", ...)
{
  concepts = get.unified.concepts(data)
  mat = matrix(0, nrow=length(concepts), ncol=length(concepts), dimnames=list(concepts, concepts))
  for (m in 1:length(data$maps))
  {
    res = get.edgelist(data$maps[[m]]$map)
    if (dim(res)[1] > 0)
      for (i in 1:dim(res)[1])
        mat[which(concepts == res[i, 1]), which(concepts == res[i, 2])] = mat[which(concepts == res[i, 1]), which(concepts == res[i, 2])] + 1
  }
  mat = mat + t(mat)
  
  keep = which(colSums(mat) > threshold)
  if (prune.edges)
  {
    mat[which(mat < threshold)] = 0
    keep = which(colSums(mat) > 0)
  }
  
  mat = mat[keep, keep]
  
  mat = data$count + 1 - mat
  res = pathfinder(mat, min(q, dim(mat)[1]-1), r)
  res[which(res == Inf)] = 0
  
  graph = graph.adjacency(res, mode="plus", weighted=T)
  if (filename != "")
  {
    f = file(filename, open="wt")
    for (i in 1:dim(res)[1])
    {
      write(paste(i, dimnames(res)[[1]][i]), f)
    }
    write("#", f)
    e = get.edgelist(graph)
    for (i in 1:dim(e)[1])
    {
      write(paste(which(dimnames(res)[[1]] == e[i, 1]), which(dimnames(res)[[1]] == e[i, 2])), f)
    }
    close(f)    
  }
  if (return.cm)
    return(conceptmap(graph))
  else
    return(graph)
}


#' Creating a Pathfinder network from an igraph object
#' 
#' \code{pathfinder} creates the Pathfinder network from a weighted graph based on \code{\link{pathfinder.matrix}}. It is a convenience method that can be called
#' on the result of a call to \code{\link{landscape}}
#' @param data An igraph object.
#' @param q The parameter q used in the Pathfinder algorithm. The resulting graph will be q-triangular.
#' @param r The parameter r used in the Pathfinder algorithm for the r-metric.
#' @param threshold A numeric value used for pruning the graph before the Pathfinder algorithm. The pruning works in conjunction with
#' the value of \code{prune.edges}.
#' @param prune.edges If TRUE, each entry of the weight matrix that is lower than \code{threshold} will be set to 0 and columns
#' with a resulting sum of 0 are removed. If FALSE, only columns of the weight matrix with a sum of less than \code{threshold} 
#' will be removed.
#' @param filename Optional. If specified, the resulting Pathfinder network will be stored in TGF format in the given file.
#' @param ... -
#' @return An igraph object that represents the Pathfinder network as a weighted graph.
#' @examples
#' #Create concept maps from three random graphs
#' require("igraph")
#' g1 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' g2 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' g3 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' 
#' #Create conceptmaps object from three conceptmap objects
#' simple_cms = conceptmaps(list(conceptmap(g1), conceptmap(g2), conceptmap(g3)))
#' pathfinder(landscape(simple_cms, result="graph", mode="undirected"))
#' @method pathfinder igraph
#'@export
pathfinder.igraph <- function(data, q=Inf, r=1, threshold=0, prune.edges=F, filename="", ...)
{
  mat = get.adjacency(as.undirected(data, mode="collapse"), type="both")  
  keep = which(colSums(mat) > threshold)
  if (prune.edges)
  {
    mat[which(mat < threshold)] = 0
    keep = which(colSums(mat) > 0)
  }
  
  mat = mat[keep, keep]
  
  res = pathfinder(as.matrix(mat), min(q, dim(mat)[1]-1), r)
  res[which(res == Inf)] = 0
  
  graph = graph.adjacency(res, mode="plus", weighted=T)
  if (filename != "")
  {
    f = file(filename, open="wt")
    for (i in 1:dim(res)[1])
    {
      write(paste(i, dimnames(res)[[1]][i]), f)
    }
    write("#", f)
    e = get.edgelist(graph)
    for (i in 1:dim(e)[1])
    {
      write(paste(which(dimnames(res)[[1]] == e[i, 1]), which(dimnames(res)[[1]] == e[i, 2])), f)
    }
    close(f)    
  }
  return(graph)
}

