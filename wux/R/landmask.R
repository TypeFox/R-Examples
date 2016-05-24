GetLandMask <- function( land.mask,
                         land.mask.name,
                         clip.data.info,
                         clip.land.mask.info ) {

  ## read land mask
  ## nc <- open.ncdf(land.mask, readunlim=TRUE)
  nc <- nc_open(land.mask, readunlim=TRUE)
  parameter.dims <- GetVariableDims(nc, land.mask.name)
  if ( length( parameter.dims) == 2 ) {
    nc.start <- c("lon"=0,"lat"=0)
    nc.count <- c("lon"=0,"lat"=0)
  } else {
    nc.start <- c("lon"=0,"lat"=0,"time"=1)
    nc.count <- c("lon"=0,"lat"=0,"time"=1)
  }

  nc.count[c("lon", "lat")] <- clip.land.mask.info$count
  nc.start[c("lon", "lat")] <- clip.land.mask.info$offset

  ## land.mask.values <-                                   ##old...ncdf
  ##   get.var.ncdf(nc, land.mask.name,
  ##                start=nc.start,
  ##                count=nc.count, verbose = FALSE)
  land.mask.values <-
    ncvar_get(nc, land.mask.name,
                 start=nc.start,
                 count=nc.count, verbose = FALSE)
  ## close.ncdf(nc)
  nc_close(nc)

  ## compare lon lat of data and land mask and rotate land mask if necessary
  lon.comp <- all.equal(round(clip.data.info$lon, digits=2),
                        round(clip.land.mask.info$lon, digits=2))
  if ( !is.logical(lon.comp) )
    lon.comp <- FALSE
  lat.comp <- all.equal(round(clip.data.info$lat, digits=2),
                        round(clip.land.mask.info$lat, digits=2))
  if ( !is.logical(lat.comp) )
    lat.comp <- FALSE
  ## check if data and land mask are on the same grid,
  ## if not, stop
  if ( !lon.comp & !lat.comp ) {
    if ( mean(clip.data.info$lon) !=
        mean(clip.land.mask.info$lon) ) {
      stop(cat("WUX ERROR in GetLandMask:", "\n",
               " DATA AND LAND MASK ARE NOT",
               "ON THE SAME GRID (LON VALUES DIFFER)", "\n"))
    }
    if ( mean(clip.data.info$lat) !=
        mean(clip.land.mask.info$lat) ) {
      stop(cat("WUX ERROR in GetLandMask:", "\n",
               " DATA AND LAND MASK ARE NOT",
               "ON THE SAME GRID (LON VALUES DIFFER)", "\n"))
    }
  }
  ## lon and lat are both the same
  if ( lon.comp & lat.comp ) {
    land.mask.values <- land.mask.values
  }
  ## lon are the same, but lat are different
  if ( lon.comp & !lat.comp ) {
    flip.lat.values <-
      c(length(clip.land.mask.info$lat[1,]):1)
    land.mask.values <- land.mask.values[,flip.lat.values]
  }
  ## lat are the same, but lon are different
  if ( !lon.comp & lat.comp ) {
    flip.lon.values <-
      c(length(clip.land.mask.info$lon[,1]):1)
    land.mask.values <- land.mask.values[flip.lon.values,]
  }
  ## lon and lat are not the same
  if ( !lon.comp & !lat.comp ) {
    flip.lat.values <-
      c(length(clip.land.mask.info$lat[1,]):1)
    flip.lon.values <-
      c(length(clip.land.mask.info$lon[,1]):1)
    land.mask.values <- land.mask.values[flip.lon.values,]
    land.mask.values <- land.mask.values[flip.lon.values,flip.lat.values]
  }

  if ( floor(max(range(land.mask.values, na.rm=T))) <= 1 )
    land.mask.values <- land.mask.values * 100
  land.mask.clip <- which(land.mask.values < 50)

  return(land.mask.clip)

}
