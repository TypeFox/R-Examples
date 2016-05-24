# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2012
# Version 1.0
# Licence GPL v3

if (!isGeneric("cellFromRowCol")) {
  setGeneric("cellFromRowCol", function(object, rownr, colnr)
    standardGeneric("cellFromRowCol"))
}

setMethod("cellFromRowCol", "RasterStackBrickTS",
          function(object, rownr, colnr) {
            if ( missing(rownr) | missing(colnr)) { stop('you must provide row and col number(s)') }
            as.vector(cellFromRowCol(object@raster,rownr=rownr,colnr=colnr))
          })
#-----------
if (!isGeneric("cellFromXY")) {
  setGeneric("cellFromXY", function(object, xy)
    standardGeneric("cellFromXY"))
}

setMethod("cellFromXY", "RasterStackBrickTS",
          function(object,xy) {
            if ( missing(xy)) { stop('you must provide XY coordinates as a vector, matrix, or SpatialPoints') }
            cellFromXY(object@raster,xy)
          })