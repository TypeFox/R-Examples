#' Determine nearest neighbors
#' 
#' \code{nn} determines the nearest neighbors for a set of observations based on the distance matrix.
#' 
#' This function can determine nearest neighbors in two ways:  1. by total count or 2. by distance.
#' If method = "c", then k determines the total number of neighbors to return for each observation.
#' If method = "d", then k determines the distance for which an observation is considered a neighbor.
#' 
#' @param d A square distance matrix for the set of coordinates of interest.
#' @param k The number of numbers to return (if \code{method = "c"}) or the distance for which observations are considered neighbors (if \code{method = "d"}).
#' @param method The method of determining the neighbors.  The default is \code{"c"}, specifying that the \code{k} nearest neighbors (the count of neighbors) for each observation should be returned.  The alternative is \code{"d"}, meaning that neighbors are determined by their distance from an observation.  In that case, two observations are neighbors if their separation distance is less or equal to \code{k}.
#' @param self A logical indicating whether an observation is a neighbor with itself.  The default is \code{FALSE}.  
#'
#' @return Returns the indexes of the nearest neighbors as a matrix if \code{method = "c"} and a list otherwise.  For each row or element of the list, the indexes are ordered from nearest to farthest.
#' @author Joshua French
#' @export
#' @examples 
#' data(grave)
#' # make distance matrix
#' d = as.matrix(dist(cbind(grave$x, grave$y)))
#' # 3 nearest neighbors
#' nnc = nn(d, k = 3, method = "c")
#' # nearest neighbors within k units of each observation
#' nnd = nn(d, k = 200, method = "d")
#' 

nn <- function(d, k, method = "c", self = FALSE)
{
  if(!is.matrix(d)) stop ("d must be a matrix of coordinates")
  if(nrow(d) != ncol(d)) stop("d should be a square matrix")
  if(!is.numeric(k) || k <= 0 || length(k) > 1) stop("k must be a single posive numeric value")
  if(!is.element(method, c("c", "d"))) stop("valid options for method are 'c' or 'd'")
  if(!is.logical(self)) stop("self must be a logica")
  nws = (self + 1) %% 2 # switch true to 0, false to 1
  if((k+nws >= ncol(d)) & method == "c") stop("k is larger than available number of neighbors")
  
  if(method == "c")
  {
    # implement slightly different behavior based on whether the observation 
    # is a neighbor to itself
    # nws was converted to numeric to shift index appropriately, depending
    # on what is wanted
    # ordered distances (column represent results for each observation)
    # select only needed number of rows
    # not that we increment the returned indexes returned by 1
    # if each observation is not a neighbor with itself
    return(t(apply(d, 1, order)[(1 + nws):(k + nws), ]))
  }else
  {
    # for each row in d: order the distances, only return the
    # ones with distances less than k.  Since distnces ordered, 
    # this will simply be the first q ordered elements,
    # where q is the number of observations in that row 
    # with distances <= k.  Adjustment made based on whether
    # whether we should include the observation itself (nws = TRUE)
    
    return(apply(d, 1,                  
                 FUN = function(x)
                 {
                   return(order(x)[(1+nws):(sum(x <= k))])
                 }))
  }
}
