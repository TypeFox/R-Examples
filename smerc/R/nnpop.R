#' Determine nearest neighbors
#' 
#' \code{nnpop} determines the nearest neighbors for a set of observations based on the distance matrix according to a population upperbound.
#' 
#' This function determines the nearest neighbors of each centroid based on the intercentroid distance.  The number of nearest neighbors is limited by the sum of the population values among the nearest neighbors.  The set of nearest neighbors can contain no more than \code{ubpop * sum(pop)} members of the population.  The nearest neighbors are ordered from nearest to farthest.  
#' 
#' @param d An \eqn{n\times n} square distance matrix containing the intercentroid distance between the \eqn{n} region centroids.
#' @param pop A vector of length \eqn{n} containing the population values of the \eqn{n} region centroids.
#' @param ubpop A proportion between 0 and 1 containing the upperbound for the proportion of total population contained collectively among a set of nearest neighbors.
#'
#' @return Returns the indexes of the nearest neighbors as a list.  For each element of the list, the indexes are ordered from nearest to farthest from each centroid.
#' @author Joshua French
#' @importFrom SpatialTools dist1
#' @export
#' @examples 
#' data(nydf)
#' d = SpatialTools::dist1(as.matrix(nydf[,c("longitude", "latitude")]))
#' nnout = nnpop(d, pop = nydf$pop, ubpop = 0.5)
#' 
nnpop = function(d, pop, ubpop)
{
  tpop = sum(pop) # total population
  # order distances for each region
  # results has each column showing order of indexes from shorter to largest distance
  od = apply(d, 2, order)
  sd = apply(d, 2, sort)
  
  # for each row of ordered distance matrix
  # sum the cumulative population size for expanding collection
  # of regions
  # return the largest collection of regions for which the cumulative
  # population is less than the desired proportion of the total popuation
  return(apply(od, 2, 
               FUN = function(x){
                 csum = cumsum(pop[x])
                 x[which(csum <= tpop*ubpop)] 
               }))
}
