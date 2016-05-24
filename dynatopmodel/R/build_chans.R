#' Construct a raster of channel locations from vector or topographic data
#'
#' @description The discretise and make.routing.table methods both require a raster defining the locations of the channel cells and the proportion of each river cell occupied by the channel.
#' A detailed river network (DRN) may be available in vector format. If not, the channel location can be inferred from a spatially-distributed metric, typically the topographic wetness index.
#'
#' @export build_chans
#' @param dem Elevation raster (DEM) using a projected coordinate system (e.g UTM) and regular grid spacing. Not required if atb raster supplied.
#' @param drn Detailed river network (DRN) in vector (ESRI Shapefile) format. Not required if atb raster supplied.
#' @param chan.width Channel width, in m, which by default will be applied across entire network.
#' @param atb Raster whose values provide a criteria for locating the channel. This is typically the value of the topographic wetness index (TWI) determined from the elevations. Should be in a projected coordinate system (e.g UTM) and regular grid spacing.
#'
#' For the TWI to be meaningful this raster should have a resolution of a least 30m. It can be calculated using the upslope.area method applied to the DEM and atb=T.
#' @param atb.thresh If drn not supplied then this specifies the threshold value above which cells are identified as containing part of the channel network
#' @param buffer If using a vector input then buffer the DRN by this width to capture all river cells.
#' @param single.chan If using a vector input then individual reach IDs are ignored and the first raster layer returned contains either 1 for a river cell or NA for a non-river cell. Otherwise the discretise method will create an entry for ecah channel ID
#' @return A two-band raster with the same dimensions as the elevation or ATB raster whose first layer comprises non-zero cells where identified with the channel and whose second layer holds the proportions of those cells occupied by the channel.
#' @references Kirkby, M. (1975). Hydrograph modelling strategies. In Peel, R., Chisholm, Michael, Haggett, Peter, & University of Bristol. Department of Geography. (Eds.). Processes in physical and human geography : Bristol essays. pp. 69-90. London: Heinemann Educational.
#' @examples
#'\dontrun{
#' # Build channel raster in two ways and compare the results (each with a
#' # nominal 2m channel width). Artificial drainage in this catchment has
#' # apparently introduced many channels with dimensions below the scale. Their
#' # existence would not have been inferred simply from examining the topography.
#' # This high-network connectivity is suggested as a cause of the unexpectedly
#' # high responsiveness of the catchment given high antecedent moisture conditions.
#'
#' require(dynatopmodel)
#' data("brompton")
#'
#' # (1) Using the wetness index
#'
#' a.atb <- upslope.area(brompton$dem, atb=TRUE)
#' chan.rast.1 <- build_chans(atb=a.atb$atb)
#'
#' # (2) using the DRN
#' chan.rast.2 <- build_chans(dem=brompton$dem, drn=brompton$drn, buff=5, chan.width=2)
#'
#' sp::plot(chan.rast.2[[1]], col="green", legend=FALSE)
#' sp::plot(chan.rast.1[[1]], col="blue", legend=FALSE, add=TRUE)
#' legend(fill=c("green", "blue"), legend=c("TWI", "DRN"), title="Method", x="bottomright")
#' }
build_chans <- function(dem=NULL,
                        drn=NULL,
                        atb=NULL,
                        chan.width=5,
                        buffer=10,
                        atb.thresh=0.8,
                        single.chan=TRUE)
{
  if(!is.null(dem) & !is.null(drn))
  {
    # build the reach multi band raster and save
    message("Building raster for channel(s)...")
    reaches <- build.reach.raster(dem, drn=gBuffer(drn, width=buffer, byid=!single.chan),
                                  chan.width=chan.width)

    reaches <- crop(reaches, dem)
    # Estimate proportion of river cells occupied by channel
  #  prop <- min(chan.width/xres(dem), 1)

  }
  else if(!is.null(atb))
  {
     # using the TWI to identify the channel
    reaches <- atb > atb.thresh*max(atb[], na.rm=T)
    # Estimate proportion of river cells occupied by channel
    prop <- min(chan.width/xres(atb), 1)
    reaches[which(reaches[]==0)] <- NA
    cellprops <- reaches*prop
    reaches <- addLayer(reaches, cellprops)
  }
  else
  {
    stop("Provide vector channel data or raster layer to locate channels")
  }

  names(reaches) <- c("chans", "chanprops")
  return(reaches)
}
