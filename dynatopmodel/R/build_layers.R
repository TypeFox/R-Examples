#' Construct basic landscape layer data for Dynamic TOPMODEL run
#'
#' @description Given an elevation raster this function will create a basic multi-band raster that can be used to run Dynamic TOPMODEL after applying a suitable discretisation.
#' It comprises the supplied elevations with the addition of upslope contributing area and topographic wetness index (TWI).
#' @export build_layers
#' @import raster
#' @param dem Elevation raster using a projected coordinate system (e.g UTM) and regular grid spacing. Should have a resolution of a least 30m for the TWI to be meaningful.
#' @param fill.sinks If TRUE (default) then run a sinkfill before calaculating the upslope area and TWI.
#' @param deg Threshold intercell slope to determine sinks (degrees).
#' @author Peter Metcalfe
#' @return A multi-band raster (stack) comprising, in order, the elevations, upslope area and topographic wetness index values.
#' @examples
#' \dontrun{
#' require(dynatopmodel)
#' data("brompton")
#'
#' # Upslope area and wetness index for Brompton catchment
#' layers <- build_layers(brompton$dem)
#'
#' sp::plot(layers, main=c("Elevation AMSL (m)", "Upslope area (log(m^2/m))", "TWI ((log(m^2/m))"))
#' }

build_layers <- function(dem,
                         fill.sinks=TRUE,
                         deg=0.1)  # sinkfill min threshold)    # default width of channel)
{
  layers <- stack(dem)

  # get the upslope area and a/tanb index. Sinks are first filled using the appropriate threshold angle
  message("Building upslope areas...")
  a.atb  <- upslope.area(dem, atb=TRUE,
                         fill.sinks=fill.sinks, deg=deg)

  layers <- addLayer(layers, a.atb)

  names(layers) <- c("dem", "a", "atb")

  return(layers)
}

