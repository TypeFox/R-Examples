#' Determine zones for flexibly shaped spatial scan test
#' 
#' \code{flex.zones} determines the unique zones to consider for the flexibly shaped spatial scan test of Tango and Takahashi (2005).
#' 
#' @param coords The number of cases in each region.
#' @param w The binary spatial adjacency matrix.
#' @param k The maximum number of regions to include in a zone.
#' @param lonlat A logical indioating whether the coordinates are longitude/latitude.  If so, the great circle distance is used in computing the nearest/neighbor distance matrix.
#' @return Returns a list of zones to consider for clustering.  Each element of the list contains a vector with the location ids of the regions in that zone.
#' @author Joshua French
#' @importFrom fields rdist.earth
#' @importFrom SpatialTools dist1
#' @importFrom smacpod nn
#' @importFrom spdep knearneigh
#' @importFrom igraph graph_from_adjacency_matrix induced_subgraph count_components
#' @importFrom parallel mclapply
#' @importFrom utils combn
#' @export
#' @references Tango, T., & Takahashi, K. (2005). A flexibly shaped spatial scan statistic for detecting clusters. International journal of health geographics, 4(1), 11.
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = cbind(nydf$longitude, nydf$latitude)
#' flex.zones(coords = coords, w = nyw, k = 3, lonlat = TRUE)
#' 
flex.zones = function(coords, w, k = 10, lonlat = FALSE)
{
  N = nrow(coords)
  
#   if(lonlat)
#   {
#     d = fields::rdist.earth(coords, coords, miles = FALSE)
#   }else
#   {
#     d = SpatialTools::dist1(as.matrix(coords))
#   }
#   
#   # determine k nearest neighbors, based on distance
#   mynn = smacpod::nn(d, k = k, method = "c", self = TRUE)},

  mynn = spdep::knearneigh(coords, k = (k - 1), longlat = lonlat)$nn
  mynn = cbind(1:N, mynn)
  
  # determine all sets of 10 elements that include the first element
  all_sets = unlist(lapply(as.list(0:(k-1)), FUN = function(x) utils::combn(2:k, x, FUN = function(y) c(1, y), simplify = FALSE)), recursive = FALSE)
  
  # for each centroid, determine which of the possible "all_sets" are connected zones
  # return the ones that are connected, return a NULL otherwise
  # microbenchmark(
  # {
  czones = unlist(parallel::mclapply(1:N, function(i)
  {
    wi = igraph::graph_from_adjacency_matrix(w[mynn[i, ], mynn[i, ]])
    zi = lapply(all_sets, function(set)
    {
      # get indices of regions in zone
      # if the number of clusters is 1, the regions are connected, we should return 
      # the indices of the regions in that zone
      if(igraph::count_components(igraph::induced_subgraph(wi, v = set)) == 1)
        return(mynn[i, ][set])
    })
    zi[!sapply(zi, is.null)]
  }), recursive = FALSE, use.names = FALSE)
  
  return(unique(parallel::mclapply(czones, sort)))
}

# mynn2 = knearneigh(coords, k = (k - 1))
# mynn2$nn = mynn[, -1]
# mynb = knn2nb(mynn2)
# microbenchmark(
#   czones = unlist(parallel::mclapply(1:N, function(i)
#   {
#     wi = igraph::graph_from_adjacency_matrix(w[mynn[i, ], mynn[i, ]])
#     zi = lapply(all_sets, function(set)
#     {
#       # get indices of regions in zone
#       # if the number of clusters is 1, the regions are connected, we should return 
#       # the indices of the regions in that zone
#       if(igraph::count_components(igraph::induced_subgraph(wi, v = set)) == 1)
#         return(sort(mynn[i, ][set]))
#     })
#     zi[!sapply(zi, is.null)]
#   }), recursive = FALSE, use.names = FALSE), times = 3L) 
### much slower than one in code
### this is implemented using functions in the spdep package
# microbenchmark(
# czones2 = unlist(parallel::mclapply(1:N, function(i)
# {
#   # wi = subset(mynb, is.element(1:N, c(i, mynn2$nn[i, ])))
#   these_sets = lapply(all_sets, function(x) c(i, mynn2$nn[i, ])[x])
#   zi = lapply(these_sets, function(set)
#   {
#     # get indices of regions in zone
#     # if the number of clusters is 1, the regions are connected, we should return 
#     # the indices of the regions in that zone
#     if(n.comp.nb(subset(mynb, is.element(1:N, set)))$nc == 1)
#     return(sort(set))
#   })
#   zi[!sapply(zi, is.null)]
# }), recursive = FALSE), times = 3L)