#' This function takes a data set and a group set and returns centroid distances
#' 
#' This is a convience function for a reduce and distance operation based upon
#'  multivariate predictor variables.
#' @param x The NxP matrix of raw data.
#' @param grouping The Nx1 vector of grouping factors.
#' @return An NxN pairwise distance matrix.
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
#' @export
centroid_distance <- function( x, grouping ){
  if( !is(grouping, "factor"))
    grouping <- factor( grouping )
  N <- dim(x)[1]
  if( length(grouping) != N )
    stop("Your grouping and predictor variables should have the same number of rows....  Come on now.")
  
  ret <- matrix(unlist(by( x, grouping, colMeans)),nrow=length(levels(grouping)), byrow=T)
  rownames( ret ) <- levels( grouping )
  return( ret )
} 