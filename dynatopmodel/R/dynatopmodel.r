#' Implementation of the Dynamic TOPMODEL hydrological model.
#'
#' @description
#' A native R implementation and enhancement of Dynamic TOPMODEL, Beven and
#' Freers (2001) extension to the semi-distributed hydrological model TOPMODEL.
#' It includes some digital terrain analysis functions for discretisation of
#' catchments by topographic indexes and other geo-referenced layers supplying
#' relevant landscape data.
#'
#' TOPMODEL (Beven & Kirkby, 1979) is a well-established and widely used
#' hydrological model that implements a spatial aggregation strategy
#' ("discretisation") in order to reduce its computational demands. Hydrological
#' similar areas identified by the discretisation procedure are referred to as
#' hydrological response units (HRUs). Beven and Freer (2001) introduced a
#' "dynamic" variant that addressed some of the limitations of the original
#' TOPMODEL but which retained its computational and parametric efficiency. In
#' particular, the original assumption of a quasi-steady water table was replaced
#' by time-dependent kinematic routing between and within HRUs.
#'
#' The new formulation allows a more flexible discretisation, variable upslope
#' drainage areas and spatially variable physical properties, allowing the
#' introduction of any type of landscape data to identify the HRUs. It retains
#' the core dynamics of the FORTRAN implementation but makes use of data storage
#' and vectorisation features of the R language to allow efficient scaling of the
#' problem domain.
#
#' The preprocessing routines supplied incorporate handling of geo-referenced
#' spatial data to allows it to integrate with modern GIS through
#' industry-standard file formats such as GEOTiff and ESRI Shapefiles.
#'
#' @seealso \code{\link{discretise}}
#' @seealso \code{\link{run.dtm}}
#'
#' @references
#' Beven, K. J. and M. J. Kirkby (1979). A physically based variable contributing area model of basin hydrology. Hydrol. Sci. Bull 24(1): 43-69.
#'
#' Beven, K. J. and J. Freer (2001). A Dynamic TOPMODEL. Hydrological Processes 15(10): 1993-2011.
#'
#' Metcalfe, P., Beven, K., & Freer, J. (2015). Dynamic TOPMODEL: A new implementation in R and its sensitivity to time and space steps. Environmental Modelling & Software, 72, 155-172.
#' @docType package
#' @name dynatopmodel
#' @import deSolve
#' @import xts
#' @import zoo
#' @import rgdal
#' @import sp
#' @import raster
#' @importFrom methods is slot slotNames
#' @importFrom grDevices colorRampPalette colours dev.cur dev.flush dev.hold dev.list dev.new dev.off dev.set jpeg rainbow
#' @importFrom graphics abline axis box grconvertY layout legend mtext par title
#' @importFrom stats end start time xtabs
#' @importFrom utils capture.output read.csv read.table setTxtProgressBar txtProgressBar write.table
#' @importFrom lubridate wday mday round_date hour yday
#' @import rgeos
NULL
#> NULL
