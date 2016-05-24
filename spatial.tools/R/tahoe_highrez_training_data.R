#' Point and polygon files for use with spatial.tools
#' @name tahoe_highrez_training
#' @docType data
#' @author Jonathan A. Greenberg \email{spatial.tools@@estarcion.net}
#' @keywords data
#' @examples
#' tahoe_highrez_training_polygons <- readOGR(
#' 	dsn=system.file("external", package="spatial.tools"),layer="tahoe_highrez_training")
#' spplot(tahoe_highrez_training_polygons,zcol="Class")
#' tahoe_highrez_training_points <- readOGR(
#' 	dsn=system.file("external", package="spatial.tools"),layer="tahoe_highrez_training_points")
#' spplot(tahoe_highrez_training_points,zcol="SPECIES")
#' tahoe_highrez_training_points_utm <- readOGR(
#' 	dsn=system.file("external", package="spatial.tools"),
#' 	layer="tahoe_highrez_training_points_utm")
#' print(projection(tahoe_highrez_training_points_utm))
#' @import rgdal
NULL