clhs.SpatialPointsDataFrame <- function(
  x, # data
  ...
  ){

  spl <- clhs(x = x@data, ...)

  if (is(spl, "cLHS_result")) {
    spl$initial_object <- x # replacing the data.frame by the SPDF object
    spl$sampled_data <- x[spl$index_samples, ]
  }
  else
    spl <- x[spl, ]
  spl
}

setMethod("clhs", "SpatialPointsDataFrame", clhs.SpatialPointsDataFrame)
