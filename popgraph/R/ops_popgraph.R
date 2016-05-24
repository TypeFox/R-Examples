#' Overload '-' operator for pairs of \code{popgraph} objects
#' 
#' An overload of the \code{-} operator for \code{popgraph} objects that
#'  removes the edges in the first one that are in the second one.
#' @param e1 A \code{popgraph} object reprenting the offspring.
#' @param e2 A \code{popgraph} object representing the parent.
#' @return A new \code{popgraph} object that represents the genotypes
#'  left over after removing the parental part (if possible).
#' @rdname popgraph-operator-minus
#' @method - popgraph
#' @export
#' @author Rodney J. Dyer \email{rjdyer@@vcu.edu}
#' @examples
#' library(igraph)
#' e1 <- as.popgraph( graph.atlas(716) )
#' e2 <- as.popgraph( graph.atlas(806) )
#' e3 <- e1 - e2
#' par(mfrow=c(1,3))
#' l <- layout.fruchterman.reingold( e1 )
#' plot(e1, layout=l)
#' plot(e2, layout=l)
#' plot(e3, layout=l)
#' par(mfrow=c(1,1))
`-.popgraph` <- function( e1, e2 ){
  if( is.na(e1) || is.na(e2))
    stop("Cannot subtract missing popgraph objects.")
  
  if( length(setdiff(V(e1)$name, V(e2)$name)) ) 
    stop("You need to have identical node sets to perform a subtraction operation.")
  
  A1 <- to_matrix( e1, mode="adjacency" )
  A2 <- to_matrix( e2, mode="adjacency" )
  
  ret <- A1 - A2
  ret[ ret < 0 ] <- 0
  ret <- as.popgraph( ret )
  return( ret )
}
