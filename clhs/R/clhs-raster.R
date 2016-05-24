clhs.Raster <- function(
  x, # data
  ...
  ){
  spdf <- rasterToPoints(x, spatial = TRUE)
  spl <- clhs(x = spdf, ...)
  
  if (is(spl, "cLHS_result"))
    spl$initial_object <- x
  else {
    if (ncol(spl) < 2)
      spl <- raster(spl)
    else
      spl <- stack(spl)
  }

  spl
}

setMethod("clhs", "Raster", clhs.Raster)
