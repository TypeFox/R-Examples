#   if(!exists("Sys.setenv", envir = baseenv())
#     Sys.setenv <- Sys.putenv
#.First.lib <- function(lib, pkg) {
##.onLoad <- function(lib, pkg) {
##  require(methods, quietly = TRUE, warn.conflicts = FALSE)
##  require(sp)
#  .rgdal_old.PROJ_LIB <- Sys.getenv("PROJ_LIB")
#  Sys.setenv("PROJ_LIB"=system.file("proj", package = "rgdal")[1])
#  .rgdal_old.GDAL_DATA <- Sys.getenv("GDAL_DATA")
#  Sys.setenv("GDAL_DATA"=system.file("gdal", package = "rgdal")[1])
#
#  library.dynam('rgdal', pkg, lib)
#
#  .Call('RGDAL_Init', PACKAGE="rgdal")
##
##
#  cat('Geospatial Data Abstraction Library ')
#  cat('extensions to R successfully loaded\n')
#  
#}
#
#.Last.lib <- function(lib, pkg) {
##.onUnload <- function(libpath) {
#    Sys.setenv("PROJ_LIB"=.rgdal_old.PROJ_LIB)
#    Sys.setenv("GDAL_DATA"=.rgdal_old.GDAL_DATA)
#}
#
