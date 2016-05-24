#' Generate a network routing table
#'
#' @description Generates a network width table for a catchment. When passed to the run.dtm routine this will be used to route channel flows to the outlet during a Dynamic TOPMODEL run.
#' @export build_routing_table
#' @details Dynamic TOPMODEL routes channel flow to the outlet by a network-width approach (see Beven, 2012, pp. 97-97). A time-delay histogram is produced using the table. When any flow is distributed to the channel "unit" it is immediately redistributed across future time steps according to the proportions found in the histogram. These flows are then added to future outputs from the model.
#' @author Peter Metcalfe
#' @param dem Elevation raster using a projected coordinate system (e.g UTM) and a regular grid spacing. Areas outside the catchment should be set to NA
#' @param breaks Number of distance intervals
#' @param chans Optional raster of the same dimensions and resolution as the DEM. Non-zero cells in this raster are considered to contain a river channel. If not supplied then flowpaths from the entire catchment area are considered.
#' @param len.fun For large rasters the flow.len function can be very slow and many paths fail to reach a single outlet cell. This applies a simple straight line distance to the outlet to obtain a rough approximation.
#' @param outlet Index of cell or cells identified with the catchment outlet
#' @return A two-column data.frame. Its first column is the average flow distance to the outlet, in m, the second the proportions of the catchment channel network within each distance category.
#'
#' @references  Beven, K. J. (2012). Rainfall-runoff modelling : the primer. Chichester, UK, Wiley-Blackwell.
#' @examples
#' \dontrun{
#' # Create a routing table for the Brompton test case and show histogram
#'
#' data(brompton)
#'
#' tab <- build_routing_table(brompton$dem,
#'   chans=brompton$reaches,
#'   breaks=5)
#' barplot(tab[,2]*100, xlab="Mean flow distance to outlet (m)",
#' ylab="Network Width %", names.arg=tab[,1])
#' }
build_routing_table<- function(dem,
                              chans=NULL,
                              outlet=NULL,       # cell or point
                              breaks=5,
                              len.fun=flow.lens)   #simple.dist.to.outlet)
{
  outlets=NULL      # (not currently used) if routing to multiple basins these are the outlet cells for each
  #	reaches <- determine.reaches(dem, drn, reaches)
  #	drn.cells <- extract.cells(dem, drn)
  #
  #	hru <- hru[[1]][]
  # default outlet is lowest part of catchment

  # extract a dem for cell containing parts of the channel
  #   if(!is.null(drn))
  #   {
  #   	drn.cells <- extract.cells(dem, gBuffer(drn, w=buff))
  #   }
  if(!is.null(chans))
  {
    # just riprian areas
    if(!is(chans, "Raster"))
    {
      stop("Raster required")
    }
    drn.cells <- which(chans[]>0)
  }
  else
  {
    # use everything
    drn.cells <- which(!is.na(dem[]))
  }

  message("Calculating flow distances...")

  # calculate flowlengths to outlet
  lens <- do.call(len.fun, list(dem=dem, src=drn.cells,  outlet=outlet))

  # just consider lengths from the riparian zone
  #lens[setdiff(1:ncell(lens), drn.cells)]<- NA
  # bin the lengths into desired number of reach lengths
  # 	if(is.null(rlens))
  # 	{
  #   	rlen. <- seq(0, max(lens[], na.rm=T), by=rlen)
  # 	}
  lens.riv <- dem -dem  + lens

  # bin the length into the given number of classes
  len.bin <- cut(lens.riv, breaks=breaks)[]

  sel <- which(!is.na(len.bin))
  # two column matrix of reach length cats and reach ids
  len.tab <- data.frame(reach=len.bin[sel], len=lens.riv[sel])
  # mean flow distances for each of the classes
  mean.lens <- sapply(split(len.tab, as.factor(len.tab$reach)), function(tab){mean(tab$len)})
  ftab <- tabulate(len.bin[sel], breaks)   #length(rlens)-1)  #, len=lens.riv[sel])
  props <- ftab/sum(ftab)

  routing <- data.frame("flow.len"=round(mean.lens), "prop"=round(props,2))
  # just the default
  if(length(outlets)==0)
  {
    return(routing)
  }
#   for(i in 1:length(outlets))
#   {
#     if(!is.null(dists.outlets))
#     {
#       # precalculated flow distance
#       lens.riv <- dists.outlets[[i]]
#     }
#     # lens <-   flow.lens.2(dem, drn.cells, outlet)
#     else
#     {
#       # build a reach table based on flow distance to outlet, using just the cells around the drn
#       # they can, however, take paths that don't follow the DRN!
#       lens <- flow.lens(dem, src=drn.cells,  outlet=outlets[i])
#       lens.riv <- dem -dem  + lens
#       lens.riv <- lens.riv[which(!is.na(lens.riv))]
#     }
#
#     len.bin <- cut(lens.riv, breaks=rlens, labels=F)
#
#     # bin the lengths into desired number of reach lengths
#     #    rlens <- seq(0, max(lens.riv[], na.rm=T), by=rlen)
#
#     #sel <- which(!is.na(len.bin))
#     # two column matrix of reach length cats and reach ids
#     len.tab <- data.frame(reach=len.bin, len=lens.riv)
#     # calculate mean flow distances for each of the classes?
#     mean.lens <- sapply(split(len.tab, as.factor(len.tab$reach)), function(tab){mean(tab$len)})
#     ftab <- tabulate(len.bin, length(rlens[-1]))  #, len=lens.riv[sel])
#     props <- ftab   / npath
#     routing <- cbind(round(routing), "prop"=round(props,2))
#
#   }
#
#   return(routing)
}
