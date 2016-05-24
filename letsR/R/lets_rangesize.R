#' Compute species' geographic range sizes 
#' 
#' @author Bruno Vilela
#' 
#' @description This function calculates species' range sizes from a 
#' PresenceAbsence object or directly from the species' shapefiles.
#' 
#' @param x A \code{\link{PresenceAbsence}} object or an \code{SpatialPolygonsDataFrame}.
#' @param species_name  Species names in the same order as in the 
#' \code{SpatialPolygonsDataFrame} (only needed if x is a \code{SpatialPolygonsDataFrame}).
#' @param coordinates "geographical" or "planar". Indicate wheter the shapefile
#' has geographical or planar coordinates(only needed if x is a 
#' \code{SpatialPolygonsDataFrame}). 
#' @param units "cell" or "squaremeter". Indicate if the size units wanted are 
#' in number of cells occupied or in square meters(only needed if x is a 
#' \code{PresenceAbsence} object).
#'
#' @return The result is a matrix with the range size of each species.
#' If the range size accounts for the earth curvature (Yes or No) or its size unit
#' may differ for each argument combination:
#' @return 1) SpatialPolygonsDataFrame & geographical = Square meters. Yes.
#' @return 2) SpatialPolygonsDataFrame & planar = same units as the coordinates. No.
#' @return 3) PresenceAbsence & cell = number of cells. No.
#' @return 4) PresenceAbsence & squaremeter = Square meters. Yes.
#' 
#' @examples \dontrun{
#' # SpatialPolygonsDataFrame & geographical
#' data(Phyllomedusa)  
#' rangesize <- lets.rangesize(x = Phyllomedusa, 
#'                             coordinates = "geographic")
#' 
#' # SpatialPolygonsDataFrame & planar
#' rangesize2 <- lets.rangesize(x = Phyllomedusa, 
#'                              coordinates = "planar")
#' 
#' # PresenceAbsence & cell
#' data(PAM)  
#' rangesize3 <- lets.rangesize(x = PAM, 
#'                              units = "cell")
#' 
#' # PresenceAbsence & squaremeter
#' rangesize4 <- lets.rangesize(x = PAM,
#'                              units = "squaremeter")
#' }
#' 
#' 
#' @export


lets.rangesize <- function(x, species_name = x@data[, 1],
                           coordinates = "geographic", 
                           units = "cell") {
  
  
  
  
  # For PresenceAbsence
  if (class(x) == "PresenceAbsence") {
    
    # Error control
    unt <- c("cell", "squaremeter")
    if (!any(units == unt)) {
      stop(paste("The units",  units,
                 "is not implemented.", 
                 "Please check the spelling."))
    }
    
    if (units == "cell") {
      Range_Size <- colSums(x$P[, -(1:2), drop = FALSE])
    }
    
    if (units == "squaremeter") {
      
      cellsize <- NULL
      ext1 <- c(-180, 180, -90, 90)
      global <- all(as.vector(extent(x[[2]])) == ext1)
      
      if (!global) {
        grid <- rasterToPolygons(x[[2]])
        cellsize <- try(areaPolygon(grid), silent=TRUE)
      }
      
      if (class(cellsize) == "try-error" | global) {
        cellsize <- values(area(x[[2]])) * 1000000
      } 
      r <- x[[2]]      
      values(r) <- cellsize      
      cellsize2 <- extract(r, x[[1]][, 1:2])
      multi <- x[[1]][, -(1:2), drop = FALSE] * cellsize2
      Range_Size <- colSums(multi)
    }
    
  }
  
  # For SpatialPolygons
  if (class(x) == "SpatialPolygonsDataFrame") {
    
    # Error control
    coords <- c("geographic", "planar")
    if (!any(coordinates == coords)) {
      stop(paste("The coordinates",  coordinates,
                 "is not implemented.", 
                 "Please check the spelling."))
    }
    
    if (coordinates == "planar") {
      Range_Size <- sapply(slot(x, "polygons"), slot, "area")
      Range_Size <- as.matrix(Range_Size)
    }
    
    if (coordinates == "geographic") {
      Range_Size <- areaPolygon(x)
      Range_Size <- as.matrix(Range_Size)
    }
    rownames(Range_Size) <- species_name
    Range_Size2 <- aggregate(Range_Size[, 1] ~ species_name,
                             FUN = sum)
    Range_Size <- Range_Size2[, 2, drop = FALSE]
    rownames(Range_Size) <- Range_Size2[, 1]
  }
  
  # Return
  Range_Size <- as.matrix(Range_Size)
  colnames(Range_Size) <- 'Range_size'
  return(Range_Size)    
}


