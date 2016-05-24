
# ----------------------------------------------------------------
# $Author: geh $
# $Date: 2012-12-06 13:40:57 +0100 (Thu, 06 Dec 2012) $
# $Rev: 243 $
# ----------------------------------------------------------------


ReadGridFile <- function(grid.filename, mask.varname = NULL, lonlat.var.name = NULL, ...) {
  ## Reads NetCDF grid file to extract coordinates and dimensions,
  ## and optionally extracts mask variable in case of a subregionsfile.
  ##
  ## Args:
  ##   grid.filename: NetCDF grid or subregion filename.
  ##   mask.varname: Mask variable name of NetCDF file when dealing with
  ##                 subregionsfile. Default is 'NULL', which means that
  ##                 no mask variable will be read in.
  ##                 Passing any value to this argument, means no lon
  ##                 lat coordinates will be read in.
  ## 
  ## Returns:
  ##   List object containing 'lon' and 'lat' matrix elements and the dimensions
  ##   'dims' when reading a gridfile.
  ##   If reading a subregionsfile 'lon', 'lat', 'dims' will be NULL, instead
  ##   the 'mask' list entry will contain the mask variable.
  ##
  ## History:
  ##   2010-10-29 | Original code (thm)
  ##   2016-01-13 | change to "ncdf4" library
  ## 
  ## TODO(thm): write a function to detect the lon and lat variable name
  ##   (maybe it is x and y)

  ## initialize
  nc.lon <- nc.lat <- nc.mask <- NULL
  ## nc.grid <- open.ncdf(grid.filename)   old ncdf
  nc.grid <- nc_open(grid.filename)
  ## read dimensions of netcdf file
  nc.dims <- GetNetcdfDims(nc.grid)
  ## lon.dims <- names(GetVariableDims(nc.grid, "lon"))
  ## names(nc.dims)[match(lon.dims, names(nc.dims))] <-
  ##   c("lon","lat")

  ## read lat long over for whole region
  ## but only if it is not a subregionfile
  if (is.null(mask.varname)) {    
    grid.filename.no.dir <-
      tail(unlist(strsplit(grid.filename, .Platform$file.sep)), n=1)
    cat("      reading constant fields:", grid.filename.no.dir, "\n")
    if ( is.null(lonlat.var.name) ) {
     ## nc.lon <-  get.var.ncdf(nc.grid, "lon")
     ##  nc.lat <-  get.var.ncdf(nc.grid, "lat")
        nc.lon <-  ncvar_get(nc.grid, "lon")
        nc.lat <-  ncvar_get(nc.grid, "lat")
    } else {
     ## nc.lon <-  get.var.ncdf(nc.grid, lonlat.var.name["longitude"])
     ##  nc.lat <-  get.var.ncdf(nc.grid, lonlat.var.name["latitude"])
        nc.lon <-  ncvar_get(nc.grid, lonlat.var.name["longitude"])
        nc.lat <-  ncvar_get(nc.grid, lonlat.var.name["latitude"])
    }
  }
  ## if we read in a subregionfile
  if (!is.null(mask.varname))
    nc.mask <-  ncvar_get(nc.grid, mask.varname)
 
   nc_close(nc.grid)

  ## in case of a rectangular grid, generate lon-lat matrix
  lon.dim <- length(dim(nc.lon))
  lat.dim <- length(dim(nc.lat))
  if (lon.dim == 1 && lat.dim == 1) {
    cat("     ", "RECTANGULAR GRID INPUT: GENERATING LON-LAT GRID\n")
    ## make grid (by using R's recycling rule)
    nc.lon <- matrix(nc.lon, ncol =  dim(nc.lat), nrow = dim(nc.lon),
                     byrow = FALSE)
    nc.lat <- matrix(nc.lat, ncol =  dim(nc.lat), nrow = dim(nc.lon),
                     byrow = TRUE)
  }
  ## to be returned
  coordinates <- list(lon  = nc.lon,
                      lat  = nc.lat,
                      mask = nc.mask,
                      dims = nc.dims)

  return(coordinates)
}

