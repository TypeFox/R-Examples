#' Basic creation of a conceptmap object
#' 
#' \code{conceptmap} creates a conceptmap object based on either an empty concept map or on (matrix) data.
#' 
#' @param x The concept map data. For NULL, an empty concept map is created, otherwise the object is cooerced into a matrix.
#' @param ... -
#' @return A conceptmap object.
#' @examples
#' empty_cm = conceptmap(NULL)
#' @method conceptmap default
#'@export 
conceptmap.default <- function(x, ...) 
{
	if (is.null(x)) 
	{
		res <- list(map=graph.empty())
		res$concepts <- c()
		res$edges <- 0
		class(res) <- "conceptmap"
		return(res)
	}	
	else
	{	
		conceptmap(as.matrix(x))
	}
}


#' Creation of a conceptmap object from matrix data
#' 
#' \code{conceptmap} creates a conceptmap object from a given matrix of a particular format (see below).
#' 
#' @param x A matrix of character vectors with at least 3 columns. Each row is of the form:
#' start, end, label, (edge attribute 1), ..., (edge attribute m).
#' Each such row will be interpreted as a directed edge from concept "start" to concept "end" with the name "label" and (optional) m additional edge attributes.
#' The colum names of \code{map.data}, if present, will be preserved as the names for the attributes.
#' @param ... -
#' @return A conceptmap object.
#' @examples
#' data = rbind(cbind("Object", "Class", "is instance of"), cbind("Class", "Attribute", "has"))
#' cm = conceptmap(data)
#' @method conceptmap matrix
#'@export 
conceptmap.matrix <- function(x, ...) 
{
	if (dim(x)[2] < 3)
		stop("Wrong dimensions. Matrix must have at least 3 columns")
  if (dim(x)[1] == 0)
    return(conceptmap.default(NULL))
		
	node.list <- union(unique(x[, 1]), unique(x[,2]))
	n <- length(node.list)
	m <- dim(x)[1]
	adj <- Matrix(0, n, n, sparse=TRUE)

	for (i in 1:m) 
		adj[which(node.list == x[i, 1]), which(node.list == x[i, 2])] <- 1
	
	map <- graph.adjacency(adj)	
	
	V(map)$name <- node.list		
	
	if (is.null(colnames(x)))
		colnames(x) <- colnames(x, do.NULL=FALSE)
	attrs <- colnames(x)[3:dim(x)[2]]
	attrs[1] <- "name"
  
	el = get.edgelist(map)
	for (i in 1:dim(el)[1]) {
	  di <- intersect(which(x[, 1] == el[i, 1]), which(x[, 2] == el[i, 2]))[1]
	  for (j in 1:length(attrs)) 
      map <- set.edge.attribute(map, attrs[j], index = i, x[di, j + 2])
	}

	res <- conceptmap(NULL)
	res$map <- map
	res$concepts <- V(map)$name
	res$edges <- length(E(map))
	return(res)
}


#' Creation of a conceptmap object from an existing graph
#' 
#' \code{conceptmap} takes an existing igraph object and tries of coerce it into a conceptmap object (encompassing the igraph object).
#' 
#' @param x An igraph object. It must have an attribute called "name" for both vertices and edges. Additional attributes are preserved for graph, vertices and edges.
#' @param strip If TRUE, nodes without adjacent edges are removed from the graph / concept map.
#' @param ... -
#' @return A conceptmap object.
#' @examples
#' #Create conceptmap from a complete graph with 5 nodes
#' require("igraph")
#' graph = graph.full(5)
#' graph = set.vertex.attribute(graph, "name", value=1:5)
#' simple_cm = conceptmap(graph)
#' @method conceptmap igraph
#'@export 
conceptmap.igraph <- function(x, strip=TRUE, ...) 
{	
	if ((!is.null(V(x)) && is.null(V(x)$name)) || (!is.null(E(x)$name) && is.null(E(x)$name)))
		stop("igraph object must have 'name' attributes for vertices and edges")
	
	if (strip)
		x <- delete.vertices(x, which(degree(x) == 0)-1)
	
	res <- conceptmap(NULL)
	res$map <- x
	res$concepts <- V(x)$name
	res$edges <- length(E(x))
	return(res)
}


#' Convert a conceptmap object to a matrix
#' 
#' \code{as.matrix} converts a conceptmap object into a matrix. The output can be fed back into \code{\link{conceptmap.matrix}}.
#' 
#' @param x A conceptmap object.
#' @param ... -
#' @return A matrix with 3 columns and one row for each proposition of the concept map.
#' The 3 columns contain the start and end node of each proposition as well as the label of the edge (as character vectors).
#' @method as.matrix conceptmap
#'@export 
as.matrix.conceptmap <- function(x, ...)
{
	res <- get.edgelist(x$map)
  if (!is.null(E(x$map)$name))
	  res = cbind(res, E(x$map)$name)
  else
    res = cbind(res, rep("", length(E(x$map))))
  return(res)
}


#' Return basic information of a conceptmap object
#' 
#' \code{summary} returns basic information about a conceptmap object
#' 
#' @param object A conceptmap object.
#' @param ... -
#' @return A list with the number of concepts, edges and components of the concept map.
#' @method summary conceptmap
#'@export 
summary.conceptmap <- function(object, ...) 
{
  res <- list(length(object$concepts), object$edges, no.clusters(object$map, mode="weak"))
  names(res) <- c("concepts", "edges", "clusters")
  return(res)
}

#' Display basic information of a conceptmap object
#' 
#' \code{print} displays basic information. For plotting, see \code{\link{plot.conceptmap}}
#' 
#' @param x A conceptmap object.
#' @param ... -
#' @return -
#' @method print conceptmap
#'@export 
print.conceptmap <- function(x, ...) 
{
  cat(paste("Concept map with", length(x$concepts), "concepts and", x$edges, "edges\n"))
}


#' Adapt list of concepts of a conceptmap object
#' 
#' \code{modify.concepts} modifes the list of concepts according to a given list.
#' This includes removing concepts and adjacent propositions as well as adding (unconnected) concepts.
#' 
#' @param x A conceptmap object.
#' @param concept.list A vector of strings that contains the list of concepts.
#' @param ... -
#' @return A conceptmap object that encompasses exactly the concepts of \code{concept.list}.
#' Concepts not originally in \code{map} are added as isolated nodes/concepts. Concepts of \code{map} that are not in \code{concept.list}
#' are removed together with their adjacent propositions.
#' @examples
#' data = rbind(cbind("Object", "Class", "is instance of"), cbind("Class", "Attribute", "has"))
#' cm = conceptmap(data)
#' modify.concepts(cm, c("Object", "Class", "Method"))
#' @method modify.concepts conceptmap
#'@export 
modify.concepts.conceptmap <- function(x, concept.list, ...)
{  
  stopifnot(class(x) == "conceptmap") 
  del <- setdiff(x$concepts, concept.list)
  mat <- as.matrix(x)
  if (length(del) > 0 && dim(mat)[1] > 0)
  {
    rows <- c()
    for (i in 1:dim(mat)[1]) 
      if ((mat[i,1] %in% del) || (mat[i,2] %in% del))
        rows <- c(rows, i)
    if (length(rows) > 0)
      mat <- mat[-rows,]
  }
  if (is.null(dim(mat)) && length(mat) > 0)
    mat = matrix(mat, nrow=1)    
  tmp <- conceptmap(NULL)
  if (!is.null(dim(mat)) && (dim(mat)[1] > 0) && !is.null(E(x$map)))
    tmp <- conceptmap(as.matrix(mat))    
  add <- setdiff(concept.list, tmp$concepts)
  ret <- conceptmap(add.vertices(tmp$map, nv = length(add), name=add), strip=FALSE)
  return(ret)
}

