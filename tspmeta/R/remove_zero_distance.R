##' Remove any duplicate cities in a tsp instance.
##'
##' @param instance [\code{tsp_instance}]\cr
##'   TSP instance object.
##' @return New TSP instance in which all duplicate cities have been
##'   removed.
##' @export 
remove_zero_distances = function(instance) {
  ## FIXME: Handle distance only instances?
  tsp_instance(coords = unique(instance$coords))
}
