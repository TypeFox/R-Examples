#' RSAGA: SAGA Geoprocessing and Terrain Analysis in R
#'
#' RSAGA provides direct access to SAGA GIS functions 
#' including, for example, a comprehensive set
#' of terrain analysis algorithms for calculating local 
#' morphometric properties (slope, aspect, curvature), 
#' hydrographic characteristics (size, height, and
#' aspect of catchment areas), and other process-related 
#' terrain attributes (potential incoming solar radiation, 
#' topographic wetness index, and more).
#' In addition, (R)SAGA provides functions for importing 
#' and exporting different grid file formats, and tools 
#' for preprocessing grids, e.g. closing gaps or filling 
#' sinks.
#' 
#' RSAGA adds a framework for creating custom-defined 
#' focal functions, e.g. specialized filter and terrain 
#' attributes such as the topographic wind shelter index, 
#' within R.  This framework can be used to apply predict
#' methods of fitted statistical models to stacks of grids 
#' representing predictor variables.  Furthermore, 
#' functions are provided for conveniently picking 
#' values at point locations from a grid using kriging 
#' or nearest neighbour interpolation.
#' 
#' RSAGA requires SAGA GIS (versions 2.0.4 - 2.2.3) 
#' are currently supported) and its user-contributed 
#' modules to be available on your computer. These 
#' can be downloaded under GPL from
#' \url{http://sourceforge.net/projects/saga-gis/}.
#' Please check the help page for 
#' \code{\link{rsaga.env}} to make sure that RSAGA 
#' can find your local installation of SAGA. You may 
#' need to 'tell' RSAGA where to find SAGA GIS.
#' 
#' Thanks to Olaf Conrad, Andre Ringeler and all the 
#' other SAGA GIS developers and contributors of 
#' this excellent geocomputing tool! Thanks to Rainer 
#' Hurling, Johan van de Wauw, Massimo Di Stefano and 
#' others for helping to adapt SAGA to and test it on
#' unix and Max OSX.
#' 
#' @references Brenning, A., 2008. Statistical geocomputing combining R and
#'  SAGA: The example of landslide susceptibility analysis with
#'  generalized additive models. In J. Boehner, T. Blaschke and
#'  L. Montanarella (eds.), SAGA - Seconds Out (= Hamburger
#'  Beitraege zur Physischen Geographie und
#'  Landschaftsoekologie, vol. 19), p. 23-32.
#'  
#'  Conrad, O., Bechtel, M., Bock, M., Dietrich, H., Fischer, E., Gerlitz, L., Wichmann,
#'  V., & Boehner, J. (2015). System for Automated Geoscientific Analyses (SAGA) v. 2.1.4.
#'  \emph{Geoscientific Model Development}, 8, 1991-2007
#'
#' @author Alexander Brenning and Donovan Bangs
#' @docType package
#' @name RSAGA-package
#' @aliases RSAGA-package
NULL
