#' Filter species' shapefiles based on its presence, origin, and season
#' 
#' @author Bruno Vilela
#' 
#' @description Filter species shapefiles by origin, presence, and seasonal type (following IUCN types: \url{http://www.iucnredlist.org/technical-documents/spatial-data}, see metadata).
#'
#' 
#' @param shapes Object of class SpatialPolygonsDataFrame (see function \code{\link{readShapePoly}} 
#' to open these files).
#' @param presence A vector with the code numbers for the presence type to be maintained.
#' @param origin A vector with the code numbers for the origin type to be maintained.
#' @param seasonal A vector with the code numbers for the seasonal type to be maintained.
#' 
#' @return The result is the shapefile(s) filtered according to the selected types. If the filters remove all polygons, the result will be NULL.
#' 
#'  
#' @details Presence codes:
#' (1) Extant, 
#' (2) Probably Extant, 
#' (3) Possibly Extant, 
#' (4) Possibly Extinct, 
#' (5) Extinct (post 1500) &
#' (6) Presence Uncertain.
#' 
#' Origin codes:
#' (1) Native, 
#' (2) Reintroduced, 
#' (3) Introduced, 
#' (4) Vagrant &
#' (5) Origin Uncertain.
#' 
#' Seasonal codes:
#' (1) Resident, 
#' (2) Breeding Season, 
#' (3) Non-breeding Season, 
#' (4) Passage &
#' (5) Seasonal Occurrence Uncertain.
#' 
#' More info in the shapefiles' metadata.
#' 
#' @seealso \code{\link{plot.PresenceAbsence}}
#' @seealso \code{\link{lets.presab}}
#' @seealso \code{\link{lets.presab.birds}}
#' 
#' 
#' @export



lets.shFilter <- function(shapes, presence = NULL, origin = NULL,
                          seasonal = NULL) {
  
  # No filtering applied control
  if (is.null(presence) & is.null(origin) & is.null(seasonal)) {
    return(shapes)
  } else {
    
    # Upper case shapes names     
    names(shapes) <- toupper(names(shapes))
    
    # Presence filter
    if (!is.null(presence)) {
      pos <- shapes$PRESENCE %in% presence
      if (length(pos) > 0) {
        shapes <- shapes[pos, ]
      } else {
        shapes <-  NULL
      }
    }
    
    # Origin filter
    if (!is.null(origin)) {
      pos2 <- shapes$ORIGIN %in% origin
      if (length(pos2) > 0) {
        shapes <- shapes[pos2, ]
      } else {
        shapes <- NULL
      }
    }
    
    # Seasonal filter
    if (!is.null(seasonal)) {
      pos3 <- shapes$SEASONAL %in% seasonal
      if (length(pos3) > 0) {
        shapes <- shapes[pos3, ]
      } else {
        shapes <- NULL
      }
    }
    return(shapes)
  }
}
