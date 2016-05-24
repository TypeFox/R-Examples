## ###########################################################
## PURPOSE: Pull information about the grid used from the
##   Models3-formatted file.  This includes information such as the
##   origin of the grid (lower left corner coordinates in grid units).
##
## INPUT:
##   file: File name of Models3-formatted file whose projection we
##     want to use.
##
## RETURNS: List containing information about the grid, including the
##   origin point of the grid (lower left coordinates in grid units),
##   projection, grid cell spacing, etc.
##
## 
## NOTES:
##   Information about grid cell size, extent of grid, etc. is stored
##   in global attributes of the Models3-formatted file.
##
##
## REVISION HISTORY:
##   Original release: Jenise Swall, 2011-06-02
## ###########################################################
get.grid.info.M3 <- function(file){

  ## Open netCDF file which has the projection we want to use..
  nc <- nc_open(file)


  ## Find out if the file is indexed according to longitude/latitude
  ## (GDTYP=1).  If so, the horizontal units are degrees ("deg").  If
  ## not, then the units are meters ("m") by default.
  grid.type <- ncatt_get(nc, varid=0, attname="GDTYP")$value
  if (grid.type==1)
    hz.units <- "deg"
  else
    hz.units <- "m"


  ## Get information about the origin (lower left coordinates in grid
  ## units).
  x.orig <- ncatt_get(nc, varid=0, attname="XORIG")$value
  y.orig <- ncatt_get(nc, varid=0, attname="YORIG")$value

  ## Get information about the horizontal grid cell size (meters).
  x.cell.width <- ncatt_get(nc, varid=0, attname="XCELL")$value
  y.cell.width <- ncatt_get(nc, varid=0, attname="YCELL")$value

  ## Number of rows and columns tells us the extent of the grid.
  ncols <- ncatt_get(nc, varid=0, attname="NCOLS")$value
  nrows <- ncatt_get(nc, varid=0, attname="NROWS")$value
  ## Get number of vertical layers.
  nlays <- ncatt_get(nc, varid=0, attname="NLAYS")$value

  ## Now form a list to hold this information about the grid.
  grid.info.list <- list(x.orig=x.orig, y.orig=y.orig,
                         x.cell.width=x.cell.width,
                         y.cell.width=y.cell.width,
                         hz.units=hz.units,
                         ncols=ncols, nrows=nrows, nlays=nlays)


  ## Close the Models3 file.
  nc <- nc_close(nc)
  rm(nc)

  ## Return the string which can be passed to project() and other
  ## functions in R package rgdal.
  return(grid.info.list)
}
