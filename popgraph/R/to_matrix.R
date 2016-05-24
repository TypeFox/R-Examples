#' Matrix conversion for Population Graph
#' 
#' This function translates a Population to a matrix representing either 
#'  the adjacency structure, the shortest path, or edge weights.
#' @param x An object of type \code{popgraph}
#' @param mode The kind of matrix to make.  At present, the following types are
#'  available:
#'    \itemize{ 
#'      \item{adjacency}{A binary matrix representing the pairs of connected nodes (default)}
#'      \item{shortest path}{The shortest path between all nodes.}
#'      \item{edge weight}{Similar to the adjacency matrix but using edge weights instead of binary
#'        values}
#'    }
#' @param ... Optional arguemnts passed on to the distance functions.
#' @return A matrix (KxK) in size (where K is the number of nodes)
#' @export
#' @author Rodney J. Dyer <rjdyer@@vcu.edu> 
to_matrix <- function( x, mode=c("adjacency","shortest path","edge weight")[1], ... ) { 
  
  if( !inherits(x,"popgraph"))
    stop(paste("This function requires a popgraph object to function. You passed a",class(x)))
  
  ret <- NULL
  if( mode=="shortest path")
    ret <- shortest.paths( x, ... )
  else if( mode=="adjacency")
    ret <- get.adjacency( x, sparse=FALSE, ...)  
  else if( mode=="edge weight")
    ret <- get.adjacency( x, attr="weight", sparse=FALSE,... ) 
  
  if( length( V(x)$name ))
    rownames(ret) <- colnames(ret) <- V(x)$name 
  
  return( ret )
}


