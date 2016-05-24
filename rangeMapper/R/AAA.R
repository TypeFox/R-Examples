.onAttach <- function(libname, pkgname) {
	dcf <- read.dcf(file=system.file("DESCRIPTION", package=pkgname) )
	packageStartupMessage(paste('This is', pkgname, dcf[, "Version"] ))
    }

# data.table, foreach, rangeMapExport 'values'
utils::globalVariables(c('i', '.', 'value', 'CON'))

#' @import sp methods RSQLite lattice ggplot2 parallel
NULL

#' @importFrom grDevices colorRampPalette
NULL

#' @importFrom stats aggregate median na.omit
NULL

#' @importFrom utils packageDescription read.csv2 read.table
NULL

#' @importFrom magrittr %>%
NULL

#' @importFrom gridExtra grid.arrange arrangeGrob
NULL

#' @importFrom foreach foreach %do% %dopar% getDoParRegistered getDoParWorkers
NULL

#' @importFrom data.table data.table setDT setattr melt set
NULL

#' @importFrom rgeos gUnionCascaded readWKT
NULL

#' @importFrom rgdal CRSargs OGRSpatialRef readGDAL readOGR  writeGDAL
NULL

#' @importFrom raster raster rasterToPolygons projection
NULL

#' @importFrom maptools getinfo.shape
NULL

#' @importFrom classInt classIntervals
NULL

#' @importFrom RColorBrewer brewer.pal.info brewer.pal
NULL

