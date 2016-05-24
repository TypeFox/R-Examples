# Copyright 2000-2001 (c) Nicholas Lewin-Koh 
# modifications 2001-2003 Roger Bivand


get.Pcent <- function(theMap) {
  .Deprecated("", package="maptools", msg="exrtact centroid from SpatialPolygons object;\nobjects other than Spatial objects defined in the sp package are deprecated")
  p.cent <- function(poly, flag) {
    cent <- .External("RshpCentrd_2d", poly, as.integer(flag),
      PACKAGE="maptools")
    cent
  }
  if (!inherits(theMap,"Map")) stop("not a Map object")
  theShapes <- theMap$Shapes

  if (attr(theShapes,'shp.type') != 'poly')
    stop("Must be a valid polygon shapelist")

  cent<-lapply(theShapes, p.cent, 0)
  return(matrix(unlist(cent), ncol=2, byrow=TRUE))
}

