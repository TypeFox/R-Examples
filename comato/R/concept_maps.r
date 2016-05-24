#' Basic creation of a conceptmaps object
#' 
#' \code{conceptmaps} creates a conceptmaps object based on either an empty set of concept map or on a list of concept maps.
#' 
#' @param x The set of concept maps. For NULL, an empty set of concept maps is created, otherwise the object is cooerced into a list.
#' @param ... -
#' @return A conceptmaps object.
#' @examples
#' empty_cms = conceptmaps(NULL)
#' @method conceptmaps default
#'@export 
conceptmaps.default <- function(x, ...) 
{
	if (is.null(x)) 
	{
		res <- list(x=list())
		res$count <- 0
		class(res) <- "conceptmaps"
		return(res)
	}	
	else
	{	
		return(conceptmaps(as.list(x)))
	}
}


#' Creation of a conceptmaps object from a list
#' 
#' \code{conceptmaps} creates a conceptmaps object from a list of conceptmap objects.
#' 
#' @param x A list of conceptmap objects.
#' @param filter If TRUE, empty concept maps (i.e. concept maps without any proposition) are not contained in the resulting set.
#' @param ... -
#' @return A conceptmaps object.
#' @examples
#' #Create concept maps from three random graphs
#' require("igraph")
#' g1 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' g2 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' g3 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' 
#' #Create conceptmaps object from three conceptmap objects
#' simple_cms = conceptmaps(list(conceptmap(g1), conceptmap(g2), conceptmap(g3)))
#' @method conceptmaps list
#'@export 
conceptmaps.list <- function(x, filter=T, ...) 
{
  filtered = c()
	for (m in x)
	{
	  if (class(m) != "conceptmap")
	    stop("List must only contain objects of class conceptmap!")
	  if (!filter || (m$edges > 0))
	    filtered = c(filtered, list(m))
	}
		
	res <- conceptmaps(NULL)
	res$maps <- filtered
	res$count <- length(filtered)
	return(res)
}


#' Creation of a conceptmaps object from a matrix
#' 
#' \code{conceptmaps} creates a conceptmaps object from a set of concept maps represented as a matrix.
#' 
#' @param x A matrix that represents a set of concept maps. The first colum is taken to identify the map, i.e.
#' for each value occuring in the first column, the rows with identical values are extracted and \code{\link{conceptmap.matrix}}
#' is called on the matrix of these rows and the remaining columns to create a conceptmap object. 
#' @param filter If TRUE, empty concept maps (i.e. concept maps without any proposition) are not contained in the resulting set.
#' @param ... -
#' @return A conceptmaps object.
#' @examples 
#' data = rbind(
#'  cbind("1", "Object", "Class", "is instance of"),
#'  cbind("1", "Object", "Attribute", "has"),
#'  cbind("2", "Class", "Attribute", "possesses"),
#'  cbind("2", "Attribute", "Data-type", "has"),
#'  cbind("3", "Object", "Class", "is instance of")
#'  )
#'  cms = conceptmaps(data)
#' @method conceptmaps matrix
#'@export 
conceptmaps.matrix <- function(x, filter=T, ...) 
{
	vals <- unique(x[,1])
	res <- list()
	for (v in vals) 
	{
		rows <- which(x[,1] == v)
		t <- x[rows, 2:dim(x)[2], drop=FALSE]
		res <- c(res, list(conceptmap(t)))
	}
	return(conceptmaps(res, filter))
}


#' Convert a conceptmaps object to a matrix
#' 
#' \code{as.matrix} converts a conceptmaps object into a matrix. The output can be fed back into \code{\link{conceptmaps.matrix}}.
#' 
#' @param x A conceptmaps object.
#' @param ... -
#' @return A matrix with 4 columns and one row for each proposition of one of the concept maps.
#' The 4 columns contain an id of the map (starting from 1) and then the start and end node of each proposition as well as the label of the edge (as character vectors).
#' @examples
#' data = rbind(
#'  cbind("1", "Object", "Class", "is instance of"),
#'  cbind("1", "Object", "Attribute", "has"),
#'  cbind("2", "Class", "Attribute", "possesses"),
#'  cbind("2", "Attribute", "Data-type", "has"),
#'  cbind("3", "Object", "Class", "is instance of")
#'  )
#'  cms = conceptmaps(data)
#'  as.matrix(cms) 
#' @method as.matrix conceptmaps
#'@export 
as.matrix.conceptmaps <- function(x, ...)
{
  res = c()
  first_col = c()
  i= 1
  for (m in x$maps)
  {
    temp = as.matrix(m)
    res = rbind(res, temp)
    first_col = c(first_col, rep(i, nrow(temp)))
  }
  return(cbind(first_col, res))
}


#' Return basic information of a conceptmaps object
#' 
#' \code{summary} returns basic information about a conceptmaps object
#' 
#' @param object A conceptmaps object.
#' @param ... -
#' @return A matrix with one column for each concept map in the set and the number of concepts, edges, and components of this map
#' respectively in 3 rows.
#' @method summary conceptmaps
#'@export 
summary.conceptmaps <- function(object, ...) 
{
  res <- matrix(0, nrow=3, ncol=object$count)
  for (m in 1:object$count) 
  {
    s <- summary(object$maps[[m]])
    res[,m] <- c(s$concepts, s$edges, s$clusters)
  }
  rownames(res) <- c("concepts", "edges", "clusters")
  colnames(res) <- colnames(res, do.NULL=FALSE, prefix="Map ")
  return(res)
}


#' Display basic information of a conceptmaps object
#' 
#' \code{print} displays basic information. For plotting, see \code{\link{plot.conceptmaps}}
#' 
#' @param x A conceptmaps object.
#' @param ... -
#' @return -
#' @method print conceptmaps
#'@export 
print.conceptmaps <- function(x, ...) 
{
  cat(paste("Collection of concept maps with", x$count, "maps\n"))
}


#' Unify sets of conceptmaps
#' 
#' \code{merge} takes two conceptsmaps objects and merges the underyling sets of conceptmaps.
#' 
#' @param x A conceptmaps object.
#' @param y A conceptmaps object.
#' @param ... -
#' @return A conceptmaps object that consist of the maps of \code{x} and \code{y}.
#' @examples
#'  data = rbind(
#'  cbind("1", "Object", "Class", "is instance of"),
#'  cbind("1", "Object", "Attribute", "has"),
#'  cbind("2", "Class", "Attribute", "possesses"),
#'  cbind("2", "Attribute", "Data-type", "has"),
#'  cbind("3", "Object", "Class", "is instance of")
#'  )
#'  cm1 = conceptmaps(data[1:2,])
#'  cm2 = conceptmaps(data[3:5,])
#'  merge(cm1, cm2)
#' @method merge conceptmaps
#'@export 
merge.conceptmaps <- function(x, y, ...)
{
  stopifnot(class(x) == "conceptmaps" && class(y) == "conceptmaps")  
  res <- conceptmaps(NULL)
  res$maps <- c(x$maps, y$maps)
  res$count <- x$count + y$count
  return(res)
}


#' Select a subset of a set of conceptmaps
#' 
#' \code{splice} selects a subset of a set of concept maps and returns them as a new conceptmaps object.
#' 
#' @param maps A conceptmaps object.
#' @param keep A numeric vector containing the indices of the maps in \code{maps} that should be retained in the subset.
#' Regular R list indexing is used.
#' @return A conceptmaps object that consist of the maps with indiced of \code{maps}.
#' @examples
#'  data = rbind(
#'  cbind("1", "Object", "Class", "is instance of"),
#'  cbind("1", "Object", "Attribute", "has"),
#'  cbind("2", "Class", "Attribute", "possesses"),
#'  cbind("2", "Attribute", "Data-type", "has"),
#'  cbind("3", "Object", "Class", "is instance of")
#'  )
#'  cms = conceptmaps(data)
#'  
#'  splice(cms, c(1,3))
#'@export
splice <- function(maps, keep) 
{
  stopifnot(class(maps) == "conceptmaps")      
  return(conceptmaps(maps$maps[keep], filter=F))
}


#' Finding all concepts used.
#' 
#' \code{get.unified.concepts} identifies the common superset of concepts that is used by the maps of a conceptmaps object.
#' 
#' @param maps A conceptmaps object.
#' @return A vector of strings of each concepts that appears in at least one of the maps of \code{maps}.
#'@export
get.unified.concepts <- function(maps) 
{
  stopifnot(class(maps) == "conceptmaps")  	
  l <- lapply(maps$maps, FUN=function(x) {return(x$concepts)})
  return(unique(unlist(l)))
}


#' Unifying the concepts of a conceptmap object.
#' 
#' \code{unify.concepts} first calls \code{\link{get.unified.concepts}} on the maps of a conceptmaps object and then calls
#' \code{\link{modify.concepts.conceptmaps}} on each of the constituent maps. Afterwards, therefore, each map of the conceptmaps object will
#' share the same common superset of concepts.
#' @param maps A conceptmaps object.
#' @return A conceptmaps object in of the same map of \code{maps}, in which every map shares the same concepts.
#' @examples
#'  data = rbind(
#'  cbind("1", "Object", "Class", "is instance of"),
#'  cbind("1", "Object", "Attribute", "has"),
#'  cbind("2", "Class", "Attribute", "possesses"),
#'  cbind("2", "Attribute", "Data-type", "has"),
#'  cbind("3", "Object", "Class", "is instance of")
#'  )
#'  cms = conceptmaps(data)
#'  
#'  unify.concepts(cms)
#'@export
unify.concepts <- function(maps)
{
  stopifnot(class(maps) == "conceptmaps")  
  l <- get.unified.concepts(maps)
  m <- lapply(maps$maps, modify.concepts, l)
  return(conceptmaps(m, filter=F))
}


#' Modifying the concepts of all maps of a conceptmaps object.
#' 
#' \code{modify.concepts} calls \code{\link{modify.concepts.conceptmap}} for each conceptmap object of a conceptmaps object.
#' Therefore, all concept maps will share the same set of concepts afterwards.
#' @param x A conceptmaps object.
#' @param concept.list A vector of strings that contains the list of concepts.
#' @param filter If TRUE, concept maps that contain no propositions after the concept modification are removed from the result.
#' @param ... -
#' @return A conceptmaps object that contains (possibly a subset of) the maps of \code{maps} in which every map contains the concepts of \code{concept.list}.
#' @examples
#'  data = rbind(
#'  cbind("1", "Object", "Class", "is instance of"),
#'  cbind("1", "Object", "Attribute", "has"),
#'  cbind("2", "Class", "Attribute", "possesses"),
#'  cbind("2", "Attribute", "Data-type", "has"),
#'  cbind("3", "Object", "Class", "is instance of")
#'  )
#'  cms = conceptmaps(data)
#'  
#'  modify.concepts(cms, c("Object", "Class"), filter=TRUE)
#' @method modify.concepts conceptmaps
#'@export
modify.concepts.conceptmaps <- function(x, concept.list, filter=F, ...)
{
  stopifnot(class(x) == "conceptmaps") 
  m <- lapply(x$maps, modify.concepts, concept.list)
  return(conceptmaps(m, filter))
}
