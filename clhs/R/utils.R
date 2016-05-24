.create_unique_colname <- function(
  id = "id", # candidate name
  nm # vector containing the existing names
  ){
  new_id <- id
  k <- 0
  while(new_id %in% nm) {
    new_id <- paste(id, k, sep = "")
    k <- k + 1
  }
  id
}

.is.spatial <- function(x) {
  if (inherits(x, "Spatial") | inherits(x, "Raster"))
    res <- TRUE
  else
    res <- FALSE

  res
}
