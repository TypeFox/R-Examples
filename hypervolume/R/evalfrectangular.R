evalfrectangular <- function(data, bandwidth, points, verbose=T)
{
  if (verbose==TRUE) {cat('Building tree...')}
  tree = kdtree_build(data)
  if (verbose==TRUE) {cat(' done.\n')}
  result = rep(NA, nrow(points))
  
  pointsfinalmax <- points + repmat(t(as.matrix(bandwidth,nrow=1)), nrow(points), 1)
  pointsfinalmin <- points - repmat(t(as.matrix(bandwidth,nrow=1)), nrow(points), 1)
  
  pointsmax_numeric = t(as.matrix(pointsfinalmax))
  pointsmin_numeric = t(as.matrix(pointsfinalmin))
  
  nr = nrow(pointsfinalmax)
  nc = ncol(pointsfinalmax)
  
  if (verbose==TRUE) {cat('Querying tree...')}
  result <- kdtree_range_query_multiple(tree, pointsmin_numeric, pointsmax_numeric, nr, nc, verbose)
  if (verbose==TRUE) {cat(' done.\n')}
  
  rm(tree); gc(reset=T); # make sure the memory is released
  
  return(result)
}