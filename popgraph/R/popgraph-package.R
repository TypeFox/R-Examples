#############################################
#        _                 _       _        #
#     __| |_   _  ___ _ __| | __ _| |__     #
#    / _` | | | |/ _ \ '__| |/ _` | '_ \    #
#   | (_| | |_| |  __/ |  | | (_| | |_) |   #
#    \__,_|\__, |\___|_|  |_|\__,_|_.__/    #
#          |___/                            #
#                                           #
#############################################

#' Package for the creation of population graph objects
#' 
#' popgraph is a pacakge that is designed to create and manipulate
#'  population graph objects (Dyer & Nason 2004, MolEcol).  This 
#'  specific package was made in conjunction with the book 
#'  "Applied Landscape Genetics" by R.J. Dyer.
#'
#' \tabular{ll}{
#' Package: \tab popgraph\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0\cr
#' Date: \tab 2013-03-06\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' There are some very good examples of the components of this package are used
#'  in the vignettes for this package.
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
#' @docType package
#' @keywords package
#' @import ggplot2 
#' @import igraph
#' @import Matrix
#' @import sampling
#' @importFrom sp Lines Line SpatialPoints SpatialLines bbox
#' @name popgraph-package
#' @rdname popgraph-package
#'
NULL

#' Sonoran desert altitude.
#' 
#' This is a raster file for altitude in the Sonoran desert
#'  retion coincident with the Lophocereus and Upiga data
#'  sets.
#' @name alt
#' @docType data
#' @keywords data
NULL

#' Metadata for Baja Populations.
#' 
#' This is metadata associated with the sampling locations
#'  for the Lophocereus and Upiga data sets.
#' @name baja
#' @docType data
#' @keywords data
NULL

#' Lophocereus population graph
#' 
#' This is the population graph for the Lophocereus data
#'  that is discussed in Dyer & Nason (2004).
#' @name lopho
#' @docType data
#' @keywords data
NULL

#' Upiga population graph
#' 
#' This is the population graph for the Upiag data
#'  that is currently unpublished.
#' @name upiga
#' @docType data
#' @keywords data
NULL


#' Araptus attenuatus data
#' 
#' This is the Araptus attenuatus data in mv format
#' @name arapat_mv
#' @docType data
#' @keywords data
NULL

#' Arapat popualtion data
#' 
#' The population strata for the arapat_mv data
#' @name arapat_pop
#' @docType data
#' @keywords data
NULL




