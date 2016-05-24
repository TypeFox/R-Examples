
# ----------------------------------------------------------------
# $Author: thm $
# $Date: 2015-04-02 15:48:58 +0200 (Thu, 02 Apr 2015) $
# $Rev: 333 $
# ----------------------------------------------------------------

plotit <- function(model.name,
                   sub.name,
                   single.sub.data.clip,
                   data.list,
                   plot.subregion) {
  ## Plots an image of the data obtained from the NetCDF file being clipped
  ## according to the current subregion. This is for diagnostic purpose only.
  ##
  ## Args:
  ##   model.name: name of model
  ##   sub.name: name of subregion
  ##   single.sub.data.clip: clipping information
  ##   data.list: list of data in subregion
  ##   plot.subregion: plotting information
  ##
  ## Returns:
  ##   plot

  ## store image as png in save.subregions.plots
  png(paste(plot.subregion$save.subregions.plots,
            model.name, "_", sub.name, ".png", sep = "" ))

  ## data of first season
  x.list <- data.list[[1]]

  ## getting xlim and ylim info
  if ( !is.null(plot.subregion$xlim) ) {
    x.lim <- plot.subregion$xlim
  } else {
    x.lim <- c(range(single.sub.data.clip$lon[!is.na(x.list)])[1] - 1,
               range(single.sub.data.clip$lon[!is.na(x.list)])[2] + 1)

  }

  if ( !is.null(plot.subregion$ylim) ) {
    y.lim <- plot.subregion$ylim
  } else {
    y.lim <- c(range(single.sub.data.clip$lat[!is.na(x.list)])[1] - 1,
               range(single.sub.data.clip$lat[!is.na(x.list)])[2] + 1)
  }

  ## get a heuristic factor, so the points are not too small
  if (!is.null(plot.subregion$cex))
    fac <- plot.subregion$cex
  ## This is experimental... 
  ## ## get a heuristic factor, so the points are not too small. We
  ## ## therefore adjust to the size of the area, relative to the amount
  ## ## of pixels.
  ## n.lats <- length(unique(single.sub.data.clip$lat[!is.na(x.list)]))
  ## n.lons <- length(unique(single.sub.data.clip$lon[!is.na(x.list)]))
  ## ## heuristic... this ratio should be 0.1 (thm made this up!)
  ## made.up.ratio <- 100
  ## lon.fac <- n.lons * (x.lim[2] - x.lim[1])/made.up.ratio
  ## lat.fac <- n.lats * (y.lim[2] - y.lim[1])/made.up.ratio
  ## ## take one of those as factor only...
  ## do.shrink <- any(c(lon.fac, lat.fac) < 1)
  ## do.magn <- any(c(lon.fac, lat.fac) > 1)
  ## if (do.shrink) {
  ##   fac <- min(c(lon.fac, lat.fac)[c(lon.fac, lat.fac) < 1])
  ## }  else if (do.magn) {
  ##   fac <- max(c(lon.fac, lat.fac)[c(lon.fac, lat.fac) > 1])
  ## }

  ## plot points according to the subregion - the circle size
  ## represents the weighting factor
  weight.cex <- single.sub.data.clip$weight[which(!is.na(x.list))]
   is.area.fraction <- "TRUE"
   if (is.null(weight.cex)){
    is.area.fraction <- "FALSE"
    weight.cex <- 1                   #if no area.fraction, then set cex to 1
  }
  plot(single.sub.data.clip$lon[!is.na(x.list)],
       single.sub.data.clip$lat[!is.na(x.list)],
       xlim = x.lim,
       ylim = y.lim,
       cex =  weight.cex * fac,         #weight times cex factor from userinput
       col = "red", pch = 20, xlab = "lon", ylab = "lat",
       main = paste("climate model:", model.name, "\nregion:", sub.name, "\narea.fraction: ",is.area.fraction,  sep =" "))

  ## plot grid and borders
  grid()
  bound.sp <- getMap(resolution="less islands")
  plot(bound.sp, add=TRUE)

  ## close device
  dev.off()
}


