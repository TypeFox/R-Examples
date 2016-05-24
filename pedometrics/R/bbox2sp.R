#' Create Spatial object from a bounding box
#' 
#' This function takes the bounding box of a Spatial* object and creates a
#' SpatialPoints* or SpatialPolygons* object from it.
#' 
#' @param obj Object of class Spatial*.
#' 
#' @param sp Class of the resulting object. Available options are
#' \code{"SpatialPoints"}, \code{"SpatialPointsDataFrame"},
#' \code{"SpatialPolygons"} and \code{"SpatialPolygonsDataFrame"}.
#' 
#' @param keep.crs Logical for assigning the same coordinate reference system
#' to the resulting Spatial* object.
#' 
#' @return An object of class SpatialPoints* or SpatialPolygons*.
#' 
#' @note Some of the solutions used to build this function were found in the
#' source code of the R-package \pkg{intamapInteractive}. As such, the authors
#' of that package, Edzer Pebesma <\email{edzer.pebesma@@uni-muenster.de}> and
#' Jon Skoien <\email{jon.skoien@@gmail.com}>, are entitled
#' \sQuote{contributors} to the R-package \pkg{pedometrics}.
#' 
#' @author Alessandro Samuel-Rosa <\email{alessandrosamuelrosa@@gmail.com}>
#' @references Edzer Pebesma, Jon Skoien with contributions from Olivier Baume,
#' A. Chorti, D.T. Hristopulos, S.J. Melles and G. Spiliopoulos (2013).
#' \emph{intamapInteractive: procedures for automated interpolation - methods
#' only to be used interactively, not included in intamap package.} R package
#' version 1.1-10.  \url{http://CRAN.R-project.org/package=intamapInteractive}
#' @keywords misc spatial
#' @export
#' @examples
#' require(sp)
#' data(meuse)
#' coordinates(meuse) <- ~ x + y
#' bbox2sp(meuse, keep.crs = FALSE)
# FUNCTION #####################################################################
bbox2sp <- 
  function (obj, sp = "SpatialPolygons", keep.crs = TRUE) {
    
    # Check if suggested packages are installed
    pkg <- c("sp")
    id <- !sapply(pkg, requireNamespace, quietly = TRUE)
    if (any(id)) {
      pkg <- paste(pkg[which(id)], collapse = " ")
      stop(paste("Package(s) needed for this function to work but not",
                 "installed: ", pkg, sep = ""), call. = FALSE)
    }
    
    if (!inherits(obj, "Spatial")) {
      stop ("obj should be of class Spatial")
    }
    
    if (keep.crs) {
      if (is.na(sp::proj4string(obj))) {
        stop ("obj DOES NOT have a coordinate reference system")
      }
    }
    
    bb <- sp::bbox(obj)
    bbx <- c(bb[1, 1], bb[1, 2], bb[1, 2], bb[1, 1], bb[1, 1])
    bby <- c(bb[2, 1], bb[2, 1], bb[2, 2], bb[2, 2], bb[2, 1])
    if (sp == "SpatialPoints") {
      bb <- sp::SpatialPoints(data.frame(bbx, bby))
      bb <- bb[1:4, ]
    }
    if (sp == "SpatialPointsDataFrame") {
      bb <- data.frame(bbx, bby)
      bb <- bb[1:4, ]
      bb <- sp::SpatialPointsDataFrame(bb, data = data.frame(ID = c(1, 2, 3, 4)))
    }
    if (sp == "SpatialPolygons") {
      bb <- sp::SpatialPoints(data.frame(bbx, bby))
      bb <- sp::Polygons(list(sp::Polygon(bb)), ID = as.character(1))
      bb <- sp::SpatialPolygons(list(bb))
    }
    if (sp == "SpatialPolygonsDataFrame") {
      bb <- sp::SpatialPoints(data.frame(bbx, bby))
      bb <- sp::Polygons(list(sp::Polygon(bb)), ID = as.character(1))
      bb <- sp::SpatialPolygons(list(bb))
      bb <- sp::SpatialPolygonsDataFrame(bb, data = data.frame(ID = 1))
    }
    if (keep.crs) {
      sp::proj4string(bb) <- sp::proj4string(obj)
    }
    return (bb)
  }
