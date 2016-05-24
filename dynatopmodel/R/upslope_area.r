require(topmodel)
require(raster)

#' Upslope contributing area and wetness index calculation
#'
#' @description Determine upslope contributing area based on an elevation raster and, optionally, compute the topographic wetness index.
#' @export upslope.area
#' @import raster
#' @import topmodel
#' @param dem   Elevation raster (in m), using a  projected coordinate system with identical x and y resolutions.
#' @param log   Return the natural log of the values.
#' @param atb   If TRUE, include both the upslope contributing area and the topographic wetness index \eqn{ln(a/tan(beta))}. Otherwise calculate just the upslope area.
#' @param deg   Minimum intercell slope to identify with a sink (degrees).
#' @param fill.sinks Fill sinks before calculation using the threshold angle given by deg.
#' @note This is a wrapper to the function implemented in the TOPMODEL package by Wouter Buytaert.
#' @author Peter Metcalfe and Wouter Buytaert
#' @references Quinn, P. F., Beven, K. J., & Lamb, R. (1995). The In (a/tan/beta) index: How to calculate it and how to use it within the Topmodel framework. Hydrological processes, 9(2), 161-182.
#' @examples
#'\dontrun{
#' require(dynatopmodel)
#' data(brompton)
#'
#' a.atb <- upslope.area(brompton$dem, atb=TRUE)
#' sp::plot(a.atb, main=c("Upslope area (log(m^2/m))", "TWI log(m^2/m)"))
#' }
upslope.area <- function(dem, log=TRUE, atb=FALSE, deg=0.1, fill.sinks=TRUE)
{
  # check
  if(xres(dem)!=yres(dem))
  {
    stop("Raster has differing x and y cell resolutions. Check that it is in a projected coordinate system (e.g. UTM) and use raster::projectRaster to reproject to one if not. Otherwise consider using raster::resample")
  }
  # any sinks still present may give strange results
  #  sink(file="e:/junk/sink.txt")
  #  on.exit(sink(NULL))
  if(fill.sinks)
  {
    # use capture.output to supress the function console output
    capture.output(dem <- invisible(raster::setValues(dem, topmodel::sinkfill(raster::as.matrix(dem), res=xres(dem), degree=deg))))
  }
  topidx <- topmodel::topidx(raster::as.matrix(dem), res=xres(dem))

  a <- raster::setValues(dem, topidx$area)
  if(log)
  {
    a <- log(a)
  }
  if(atb)
  {
    atb <- raster::setValues(dem, topidx$atb)
    # add the topographic index ln(a/tanB)
    a <- addLayer(a, atb)
    names(a)<-c("a", "atb")
  }
  return(a)
}

# use function from TOPMODEL
flow.lens <- function(dem,
                      src=NULL,  #  starting cells, defaults to all cells i dem
                      agg=1, # initial aggregation factor
                      max.agg=4,

                      outlet=NULL)  # A vector containing the row and column indices of the pixel representing the catchment outlet. if a single value, treated as the index of a DEM cell
{
  lens <- raster::setValues(dem, NA)
  if(length(outlet)>0)
  {
    outlet.sp <- xyFromCell(dem, outlet, spatial=T)
  }

  dem.agg <- dem
  while(agg <= max.agg & max(c(0,lens[]), na.rm=T)==0)
  {
    if(agg>1)
    {
      message("Trying a aggregated dem to determine flow lengths")
      # try a coarser
      dem.agg <- raster::aggregate(dem, agg)
      #	reaches <- aggregate(reaches, )
    }

    if(length(outlet)>0)
    {
      outlet <- extract.cells(dem.agg, outlet.sp)
      iout<- rowColFromCell(dem.agg, outlet)
    }
    else{iout <- NA}

    dem.agg <- fill.sinks(dem.agg, deg=0.1)
    lens <- raster::setValues(dem.agg, flowlength(as.matrix(dem.agg), outlet=iout))
    if(!is.null(src))
    {
      lens[setdiff(1:ncell(dem), src)]<- NA
    }
    agg <- agg+1
  }
  agg <- agg-1
  # disaggregate
  if(agg>1)
  {
    message("Disaggregating back to original resolution...")
    lens <- raster::disaggregate(lens, agg, method="bilinear")
  }
  return(raster::xres(dem)*lens)
}

extract.cells <- function(dem, drn,...)
{
  target <- extract(dem, drn, cellnumbers=T,...)
  if(is.list(target))
  {
    return(do.call(rbind, target)[,1])
  }
  else
  {
    return(target[,1])
  }
}


fill.sinks <- function(dem, deg=0.01,
         silent=T,
         ipass=1,   # perform sinkfill a maximum of this times or until all sinks filled
         fail.if.not.complete=F)
{
  DEM <- as.matrix(dem)
  res <- xres(dem)
 # stopifnot(is(DEM, "matrix"))
 # if (min(as.vector(DEM[!is.na(DEM)])) < -9000)
#    stop("DEM contains unrealistic values (< -9000)")
  #   DEM[is.na(DEM)] <- -9999
#  nrow <- dim(DEM)[1]
#  ncol <- dim(DEM)[2]

  i <- 1
  sinks.remain <- T
  while(i <= ipass & sinks.remain)
  {
    prev <- DEM
  #  DEM[is.na(DEM)] <- -9999
    capture.output(DEM <- topmodel::sinkfill(DEM, res, deg))
    diff <- sum(DEM[]-prev[], na.rm=T)
    if(diff==0)
    {
      # there are definitely no sinks left now
      sinks.remain <- F
    }
#
#       .C("sinkfill", PACKAGE = "topmodel", as.double(DEM),
#                  result = double(nrow * ncol + 2), as.integer(nrow), as.integer(ncol),
#                  as.double(res), as.double(deg))$result
#    result[result > 999998] <- NA
#     DEM <- matrix(result[3:(nrow * ncol + 2)], nrow = nrow)
#     # 100 is max number of iterations, so if reached thsi then have to run the sinkfill again
#     if(result[1]< 100 & result[1]>0)
#     {
#       sinks.remain <- F
#
#     }
#     else if(!silent)
#     {
#       cat("Sinkfill pass #", i, " No. of sinks remaining = ", result[2], "\n")
#
#     }

    i <- i + 1
  }
#   if (result[1] == -1)
#   {
#     warning("incomplete sink removal")
#   }
#   if (result[1] == 100)
#   {
#     msg <- paste("Maximum no. iterations reached (100). No. sinks remaining=", result[2])
#     message(msg)
#   }

  #  mat <- matrix(result[3:(nrow * ncol + 2)], nrow = nrow)

  return(setValues(dem, DEM))
}
