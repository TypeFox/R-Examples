#' Lake Morphometry in R
#' 
#' Lakemorpho provides a number of functions to calculate a standard suite of 
#' lake morphometry metrics.  Most of the metrics are measurements of the 
#' shape of the lake.  Metrics that rely on depth have traditionally been 
#' calculated with bathymetry data.  In the absence of bathymetry data it is
#' possible to estimate maximum depth from surrounding topography.  Lakemorpho
#' uses this approach to also estimate maximum depth, mean depth, and volume.
#' 
#' This development version of this package is available at 
#' https://github.com/USEPA/lakemorpho 
#' 
#' @references Florida LAKEWATCH (2001). A Beginner's guide to water management
#'             - Lake Morphometry (2nd ed.). Gainesville: Florida LAKEWATCH, 
#'             Department of Fisheries and Aquatic Sciences.
#'             \href{http://edis.ifas.ufl.edu/pdffiles/FA/FA08100.pdf}{Link}
#' @references Hollister, J. W., W.B. Milstead (2010). Using GIS to Estimate 
#'             Lake Volume from Limited Data. Lake and Reservoir Management. 
#'             26(3)194-199.
#'             \href{http://dx.doi.org/10.1080/07438141.2010.504321}{Link}
#' @references Hollister, J. W., W.B. Milstead, M.A. Urrutia (2011). Predicting 
#'             Maximum Lake Depth from Surrounding Topography. PLoS ONE 6(9).
#'             \href{http://dx.doi.org/10.1371/journal.pone.0025764}{link}           
#'             
#' @name lakemorpho
NULL 
