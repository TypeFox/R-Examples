#' Crop a PresenceAbsence object based on an input shapefile
#' 
#' @author Bruno Vilela
#' 
#' @description Crop a PresenceAbsence object based on a shapefile provided by the user.
#' 
#' @param x A \code{\link{PresenceAbsence}} object.
#' @param shp Object of class SpatialPolygonsDataFrame (see function \code{\link{readShapePoly}}) to crop the PresenceAbsence object.
#' @param remove.sp Logical, if \code{TRUE} the final matrix will not contain species that do not match any cell in the grid.
#' 
#' 
#' @return The result is an object of class PresenceAbsence croped.
#' 
#' 
#' @seealso \code{\link{plot.PresenceAbsence}}
#' @seealso \code{\link{lets.presab.birds}} 
#' 
#' @examples \dontrun{
#' data(PAM)
#' # PAM before crop
#' plot(PAM, xlab = "Longitude", ylab = "Latitude",
#'      main = "Phyllomedusa species richness")
#' 
#' # Crop PAM to Brazil
#' require(maptools)    
#' data(wrld_simpl)  # World map
#' Brazil <- wrld_simpl[wrld_simpl$NAME == "Brazil", ]  # Brazil (polygon)
#' PAM_crop <- lets.pamcrop(PAM, Brazil, remove.sp = TRUE)
#' plot(PAM_crop, xlab = "Longitude", ylab = "Latitude",
#'      main = "Phyllomedusa species richness (Brazil crop)",
#'      col = colorRampPalette(c("darkgreen", "yellow", "blue")))
#' }
#' 
#' @export


lets.pamcrop <- function(x, shp, remove.sp = TRUE) {
  
  # Set NA to raster data outside shp
  remover1 <- extract(x[[2]], shp, cellnumbers = T, 
                      weights = T, small = T)
  remover2 <- do.call(rbind.data.frame, remover1)[, 1]  
  values(x[[2]])[-remover2] <- NA
  
  if (all(is.na(values(x[[2]])))) {
    stop("No species left after croping the PAM")
  }
  
  # Remove cells from the matrix
  manter <- extract(x[[2]], x[[1]][, 1:2])
  x[[1]] <- x[[1]][!is.na(manter), ]
  
  # Remove species without presence
  if (remove.sp) {
    x[[1]] <- .removeSp(x[[1]])
  }
  
  # Rename columns
  x[[3]] <- colnames(x[[1]])[-(1:2)]
  
  return(x)
}
