setMethod("cellFromRowCol", signature(object = "Speclib"), 
          definition = function(object, rownr, colnr)
{
  if (object@spectra@fromRaster)
  {
    return(cellFromRowCol(object@spectra@spectra_ra, rownr, colnr))
  } else {
    stop("Speclib does not contain spectra from *raster-object")
  }
}
)

setMethod("cellFromRowColCombine", signature(object = "Speclib"), 
          definition = function(object, rownr, colnr)
{
  if (object@spectra@fromRaster)
  {
    return(cellFromRowColCombine(object@spectra@spectra_ra, rownr, colnr))
  } else {
    stop("Speclib does not contain spectra from *raster-object")
  }
}
)

setMethod("cellFromRow", signature(object = "Speclib"), 
          definition = function(object, rownr)
{
  if (object@spectra@fromRaster)
  {
    return(cellFromRow(object@spectra@spectra_ra, rownr))
  } else {
    stop("Speclib does not contain spectra from *raster-object")
  }
}
)

setMethod("cellFromCol", signature(object = "Speclib"), 
          definition = function(object, colnr)
{
  if (object@spectra@fromRaster)
  {
    return(cellFromCol(object@spectra@spectra_ra, colnr))
  } else {
    stop("Speclib does not contain spectra from *raster-object")
  }
}
)

setMethod("colFromX", signature(object = "Speclib"), 
          definition = function(object, x)
{
  if (object@spectra@fromRaster)
  {
    return(colFromX(object@spectra@spectra_ra, x))
  } else {
    stop("Speclib does not contain spectra from *raster-object")
  }
}
)

setMethod("rowFromY", signature(object = "Speclib"), 
          definition = function(object, y)
{
  if (object@spectra@fromRaster)
  {
    return(rowFromY(object@spectra@spectra_ra, y))
  } else {
    stop("Speclib does not contain spectra from *raster-object")
  }
}
)

setMethod("cellFromXY", signature(object = "Speclib"), 
          definition = function(object, xy)
{
  if (object@spectra@fromRaster)
  {
    return(cellFromXY(object@spectra@spectra_ra, xy))
  } else {
    stop("Speclib does not contain spectra from *raster-object")
  }
}
)

setMethod("cellFromLine", signature(object = "Speclib"), 
          definition = function(object, lns)
{
  if (object@spectra@fromRaster)
  {
    return(cellFromLine(object@spectra@spectra_ra, lns))
  } else {
    stop("Speclib does not contain spectra from *raster-object")
  }
}
)

setMethod("cellFromPolygon", signature(object = "Speclib"), 
          definition = function(object, p, weights=FALSE)
{
  if (object@spectra@fromRaster)
  {
    return(cellFromPolygon(object@spectra@spectra_ra, p, weights = weights))
  } else {
    stop("Speclib does not contain spectra from *raster-object")
  }
}
)

setMethod("fourCellsFromXY", signature(object = "Speclib"), 
          definition = function(object, xy, duplicates=TRUE)
{
  if (object@spectra@fromRaster)
  {
    return(fourCellsFromXY(object@spectra@spectra_ra, xy, duplicates = duplicates))
  } else {
    stop("Speclib does not contain spectra from *raster-object")
  }
}
)

setMethod("readAll", signature(object = "Speclib"), 
          definition = function(object)
{
  if (object@spectra@fromRaster)
  {
    spectra(object) <- readAll(object@spectra@spectra_ra)
  }   
  return(object)
}
)