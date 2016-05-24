##' Gets Pairs from Distances
##' 
##' For each location, we wish to determine the k closest locations.  This
##' function takes the distance matrix and computes the returns a list of the
##' k closest locations to each individual location.
##' 
##' @param dist The distance matrix describing the distance between locations.
##' @param k The number of closest neighbors to be located for each lcoation.
##' 
##' @return A named list.  Each element of the list corresponds to a particular
##' location, and the value at element i is the k closest locations to location
##' i.
##' 

getPairs = function(dist, k){
  pairs = lapply(1:nrow(dist), function(i){
    x = dist[i,]
    filt = rank(x)>1 & rank(x)<=k+1.5
    colnames(dist)[filt]
  })
  names(pairs) = colnames(dist)    
  return(pairs)
}