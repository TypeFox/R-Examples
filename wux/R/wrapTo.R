
# ----------------------------------------------------------------
# $Author: geh $
# $Date: 2012-12-06 13:40:57 +0100 (Thu, 06 Dec 2012) $
# $Rev: 243 $
# ----------------------------------------------------------------

wrapTo180 <- function(poly.sp) {
  ## Wraps all longitude values to -180 to 180 degrees.
  ##
  ## Args:
  ##   poly.sp:
  ##
  ## Returns:
  ##   object of same class as 'poly.sp', but with lon coordniates ranging
  ##   from -180 to 180 degrees.
  ##
  ## History:
  ##   2011-02-16 | Original code (thm)
  ##   2011-05-02 | major bugfix for class(poly.sp) and matrix class
  ##                handling added (thm)
  
  if (!class(poly.sp) %in% c("SpatialPolygons", "numeric", "matrix")) 
    stop("CLASS OBJECT ", class(poly.sp), " UNKNOWN TO wrapTo180, BUT FEEL FREE TO EXPAND THIS FUNCTION :)")
  
  if (class(poly.sp) == "SpatialPolygons") {
    poly.sp2 <- poly.sp
    coord <- lapply(poly.sp2@polygons[[1]]@ Polygons, function(x) {
      slot(x, "coords")})
    coord <- lapply(coord, function(x) {
      if (any(x[,1] > 180)) {
        x[,1][x[,1] > 180] <-  x[,1][x[,1] > 180] - 360}
      else {
        x <- x}
      return(x)})
    single.id <- poly.sp@polygons[[1]]@ID
    list.Polygon <- lapply(coord, Polygon) #coerce coordinates to list of Polygon objects
    poly.sp <- SpatialPolygons(list(Polygons(list.Polygon, ID = single.id)))
  }
  
  if (class(poly.sp) %in% c("numeric", "matrix")) {
    indices <- poly.sp > 180
    poly.sp[indices] <- poly.sp[indices] - 360
  }

  return(poly.sp)
}

wrapTo360 <- function(poly.sp) {
  ## Wraps all longitude values to 0 to 360 degrees.
  ##
  ## Args:
  ##   poly.sp:
  ##
  ## Returns:
  ##   object of same class as 'poly.sp', but with lon coordniates ranging
  ##   from 0 to 360 degrees.
  ##
  ## History:
  ##   2011-02-16 | Original code (thm)
  ##   2011-05-02 | major bugfix for class(poly.sp) and matrix class
  ##                handling added (thm)
  
  if (!class(poly.sp) %in% c("SpatialPolygons", "numeric", "matrix")) 
    stop("CLASS OBJECT ", class(poly.sp), " UNKNOWN TO wrapto360, BUT FEEL FREE TO EXPAND THIS FUNCTION :)")
  
  if (class(poly.sp) == "SpatialPolygons") {
    poly.sp2 <- poly.sp
    coord <- lapply(poly.sp2@polygons[[1]]@ Polygons, function(x) {
      slot(x, "coords")})
    coord <- lapply(coord, function(x) {
      if (any(x[,1] < 0)) {
        x[,1][x[,1] < 0] <-  x[,1][x[,1] < 0] + 360}
      else {
        x <- x}
      return(x)})
    single.id <- poly.sp@polygons[[1]]@ID
    list.Polygon <- lapply(coord, Polygon) #coerce coordinates to list of Polygon objects
    poly.sp <- SpatialPolygons(list(Polygons(list.Polygon, ID = single.id)))
  }
  
  if (class(poly.sp) %in% c("numeric", "matrix"))
    poly.sp[poly.sp < 0] <- poly.sp[poly.sp < 0] + 360
  
  return(poly.sp)
}
