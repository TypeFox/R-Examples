
## ----------------------------------------------------------------
## $Author: geh $
## $Date: 2015-03-31 14:00:14 +0200 (Tue, 31 Mar 2015) $
## $Rev: 331 $
## ----------------------------------------------------------------


GetAggregatedSpatioTemporalData <- function(filenames, model.name = NULL,
                                            plot.subregion,
                                            grid.filenames = NULL,
                                            subregions = NULL, parameter.name,
                                            use.land.mask, land.mask,
                                            land.mask.name,
                                            interpolate.to.grid = FALSE,
                                            temporal.aggr, startdate, enddate,
                                            ...) {
  ## Reads in all data and aggregates them according to first and second
  ## aggregation statistic and optionally interpolates to common grid.
  ##
  ## Args:
  ##   filenames: Character vector containing all filenames to be read in.
  ##   grid.filenames: Lat- lon gridfile
  ##   subregions: List of subregions to be extracted.
  ##   parameter.name: Short parameter name (with longname in "names" attribute)
  ##   interpolate.to.grid: Grid where to interpolate to.
  ##   first.aggr: List of first statistic options from user.input.
  ##   second.aggr:  List of second statistic options from user.input.
  ##   '...': Keywords for ReadNetCdfTimeData, ReadNetCdfData and
  ##          AggregateFirstStatistic
  ##
  ## Returns:
  ##   List of 2-dimensional arrays stored according to second statistic
  ##   and subregions.
  ##
  ## History:
  ##   2010-10-25 | Original code (thm)

  if (interpolate.to.grid == FALSE) {
    ## no interpolation
    aggr.data.list <-
      GetAggregatedSpatioTemporalData.nointerpol(filenames,  model.name,
                                                 plot.subregion,
                                                 grid.filenames, subregions,
                                                 parameter.name, use.land.mask,
                                                 land.mask, land.mask.name,
                                                 temporal.aggr, startdate,
                                                 enddate, ...)
  } else {
    ## no interpolation
    aggr.data.list <-
      GetAggregatedSpatioTemporalData.interpol(filenames,  model.name,
                                               plot.subregion,
                                               grid.filenames, subregions,
                                               parameter.name,
                                               interpolate.to.grid = FALSE,
                                               temporal.aggr, startdate,
                                               enddate, ...)
  }

  return(aggr.data.list)
}



GetAggregatedSpatioTemporalData.nointerpol <- function(filenames,
                                                       model.name = NULL,
                                                       plot.subregion,
                                                       grid.filenames = NULL,
                                                       subregions = NULL,
                                                       parameter.name,
                                                       use.land.mask,
                                                       land.mask,
                                                       land.mask.name,
                                                       temporal.aggr,
                                                       startdate,
                                                       enddate,
                                                       lonlat.var.name,
                                                       what.timesteps,
                                                       ...) {

  ## Reads in all data and aggregates them according to first and second
  ## aggregation statistic.
  ##
  ## Args:
  ##   filenames: Character vector containing all filenames to be read in.
  ##   grid.filenames: Lat- lon gridfile
  ##   subregions: List of subregions to be extracted.
  ##   parameter.name: Short parameter name (with longname in "names" attribute)
  ##   first.aggr: List of first statistic options from user.input.
  ##   second.aggr:  List of second statistic options from user.input.
  ##   does.plot.subregions: Boolean. If TRUE, interactive plots will be
  ##                         displayed showing a rough pixel map with the clipped
  ##                         areal data of the NetCDF file.
  ##   '...': Keywords for ReadNetCdfTimeData, ReadNetCdf and
  ##          AggregateFirstStatistic
  ##
  ## Returns:
  ##   List of 2-dimensional arrays stored according to second statistic
  ##   and subregions.
  ##
  ## History:
  ##   2010-10-25 | Original code (thm)
  ##   2010-11-24 | Read all subregion information before reading in data (thm)
  ##   2011-09-22 | added case filenames = NA (thus not reading in any data)
  ##   2011-11-30 | n.pixels calculated correctly (at least I hope)
  ##   2012-12-05 | added advanced plotting of the subregion (geh)
  ##   2012-12-05 | added land mask feature (geh)

### GET SUBREGION ATTRIBUTES
  ## get offset, count and mask on model grid (grid.filenames)
  ## for all subregions
  cat("    GETTING INFORMATION FOR ALL SUBREGIONS", "\n", sep = "")
  clip.data.info <- GetSubregionShapes(subregions = subregions,
                                       grid.filenames = grid.filenames,
                                       lonlat.var.name = lonlat.var.name,
                                       ...)
  ## get offset, count and mask for the land mask
  if ( use.land.mask ) {
    clip.land.mask.info <- GetSubregionShapes(subregions = subregions,
                                              grid.filenames = land.mask,
                                              lonlat.var.name = lonlat.var.name,
                                              ...)
  }

  ## generate empty list containing corresponding subregions for data
  if (!is.null(subregions)) {
    aggr.data.list <- vector("list", length(names(clip.data.info)))
    names(aggr.data.list) <- names(clip.data.info)
    ## warning...
    if (length(names(clip.data.info)) == 0)
      warning("YOU SEEM DONT SEEM TO HAVE GIVEN YOUR SUBREGIONS PROPER NAMES, YOU SHOULD CHANGE THAT IN YOUR user.input FILE.")
  }

  ## GET TIME ATTRIBUTES
  ## get time vector and time attributes in case filenames have been passed
  if ( !all(is.na(filenames)) ) {
    nc.time <- ReadNetCdfTimeData(filenames, what.timesteps=what.timesteps, ...)
    ## getting offset and count
    nc.time <- GetTimeCountsOffsets(nc.time, startdate, enddate, check.period = NULL)
  } else {
    nc.time <- NA
  }

### GET CLIMATE DATA AND PROCESS THEM
  ## read in and aggregate data for every subregion
  for (ii in seq(along = clip.data.info)) {

    ## extract first subregion
    single.sub.data.clip <- clip.data.info[[ii]]
    if ( use.land.mask ) {
      single.sub.land.mask.clip <- clip.land.mask.info[[ii]]
    }

    ## get name of subregion
    sub.name <- names(clip.data.info)[[ii]]

    ## read in netCDF file and return 3-dim array
    if ( !all(is.na(filenames)) ) {
      cat("\n    READ DATA FIELDS FOR SUBREGION \"", sub.name, "\"", sep = "")
      data.array <- ReadNetCdf(filenames, parameter.name,
                               time.vectors = nc.time,
                               offset = single.sub.data.clip$offset,
                               count =  single.sub.data.clip$count,
                               indices.exclude =
                               single.sub.data.clip$mask.indices.outside.subregion,
                               ...)

    } else {
      data.array <- NA
    }

    ## prompt number of pixels used
    if ( !all(is.na(filenames)) ) {
      pixels.used <- single.sub.data.clip$count.lon.lat[1] *
        single.sub.data.clip$count.lon.lat[2] -
          length(single.sub.data.clip$mask.indices.outside.subregion)
      cat("\n      ", pixels.used,
          " pixels used in this subregion", "\n", sep = "")
      cat("\n    DATA AGGREGATION FOR SUBREGION \"", sub.name, "\"", sep = "")
      cat("\n      array size before temporal aggregation: ",
          object.size(data.array) / 1024 / 1024 , " MB", "\n", sep="")
    } else {
      message("      no data read in")
    }

    ## get land mask, return values to be clipped and clip the data
    if ( use.land.mask ) {
      land.mask.clip <- GetLandMask(land.mask = land.mask,
                                    land.mask.name = land.mask.name,
                                    clip.data.info = clip.data.info[[ii]],
                                    clip.land.mask.info = clip.land.mask.info[[ii]])

      data.array <- array(apply(data.array, 3,
                                function(x){is.na(x) <- land.mask.clip; x}),
                          dim = dim(data.array))
    }

    ## temporal aggregation
    data.list <- AggregateTemporal(data.array, temporal.aggr, nc.time,
                                   startdate, enddate, what.timesteps, ...)

    ## do diagnostic plotting of the grid points in a certain area
    if ( !is.null(plot.subregion$save.subregions.plots) &
        !all(is.na(filenames)) ) {
      plotit(model.name, sub.name, single.sub.data.clip,
             data.list, plot.subregion)
    }

    cat("     array size after temporal aggregation: ",
        object.size(data.list) / 1024 / 1024 , " MB", "\n", sep="")

    aggr.data.list[[ii]]$data <- data.list
    aggr.data.list[[ii]]$weight <- single.sub.data.clip$weight
    rm(data.list)

    gc()
  }

  return(aggr.data.list)
}



GetAggregatedSpatioTemporalData.interpol <- function(filenames,
                                                     grid.filenames = NULL,
                                                     subregions = NULL,
                                                     parameter.name,
                                                     first.aggr,
                                                     second.aggr, ...)  {
  ## Reads in all data and aggregates them according to first and second
  ## aggregation statistic and interpolates to common grid.

  stop("INTERPOLATION BRANCH NOT IMPLEMENTED YET")
}



AggregateWeightedSubregions <- function( x, weight=NULL ) {
  ## helper function for aggregating spatially with weights
  ##
  ## Args:
  ##   x: list of cell values (in end-effect)
  ##   weights: their weights
  ##
  ## Returns:
  ##   list of aggregated weighted subregion values
  ##
  ## History:
  ##   2011-07-06 | Original code (mir & thm)

  if ( length(which(!is.na(x))) > 0 ) {
    if (is.null(weight)) {
      return(mean(x, na.rm=T))
    } else {
      x.na <- is.na(x)
      weight[x.na] <- NA
      weight.sum = sum(weight, na.rm=T)
      if ( weight.sum == 0 ) {
        return (NA)
      } else {
        return (sum(x*weight, na.rm=T)/weight.sum)
      }
    }
  } else {
    return(NA)
  }

}
