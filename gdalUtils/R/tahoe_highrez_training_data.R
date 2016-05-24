#' Point and polygon files for use with gdalUtils
#' @name tahoe_highrez_training
#' @docType data
#' @author Jonathan A. Greenberg \email{gdalUtils@@estarcion.net}
#' @keywords data
#' @examples
#' \dontrun{
#' tahoe_highrez_training_polygons <- readOGR(
#' 	dsn=system.file("external", package="gdalUtils"),layer="tahoe_highrez_training")
#' spplot(tahoe_highrez_training_polygons,zcol="Class")
#' tahoe_highrez_training_points <- readOGR(
#' 	dsn=system.file("external", package="gdalUtils"),layer="tahoe_highrez_training_points")
#' spplot(tahoe_highrez_training_points,zcol="SPECIES")
#' }
NULL