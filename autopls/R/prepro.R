prepro <- function (X, method = 'bn')
{
  
  method <- match.arg (method, c('bn', 'msc'), several.ok = TRUE) # more2come ..
  
  if ('bn' %in% method)
  {
    if (is.vector (X)) X <- X / sqrt (sum (X ^ 2))
    if (is.matrix (X)) X <- X / sqrt (rowSums (X ^ 2))         
    if (class (X) == 'RasterBrick' || class (X) == 'RasterStack')
      X <- X / sqrt (raster::stackApply (X ^ 2, rep (1, raster::nlayers (X)), sum))
  }
  invisible (X)
}
