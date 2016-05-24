#' Landslide Inventory and DEM
#'
#' Landslide initiation points in the \emph{Reserva Biologica San Francisco} (RBSF) area of the tropical 
#' Andes in Ecuador. The landslide inventory was mapped by Stoyan (2000) in the field and by the presence
#' of landslide scars in aerial imagery. The 10 m x 10 m digital elevation model (DEM) was triangulated
#' from aerial imagery as described by Jordan \emph{et al.} (2005) and provided courtesy of Lars Ungerechts (2010).
#' 
#' Loading this dataset also loads the object \code{dem}. Existing objects named \code{dem} may be overwritten.
#' 
#' Landslide data provided here are a subset of that used to build generalized additive models (GAMs) as landslide
#' susceptibility models by Muenchow \emph{et al.} (2012). These data correspond to those in the "natural" part of
#' the \emph{RBSF} area. Please refer to the accompanying vignette for an introductory tutorial on the use of the
#' RSAGA package for terrain analysis, geoprocessing, and model-building using these data.
#' @name landslides
#' @aliases dem
#' @format
#' A data frame of 1535 rows and 3 variables:
#' \itemize{
#' \item{\code{x}}: the x coordinate of the sample point
#' \item{\code{y}}: the y coordinate of the sample point
#' \item{\code{lslpts}}: \code{TRUE} or \code{FALSE} of landslide observation at sample point
#' }
#' 
#' A digital elevation model as a .Rd grid \code{dem}
#' \itemize{
#' \item \code{$header} list of 8 DEM header properties
#' \item \code{$data} grid elevation values (m ASL)
#' }
#' 
#' @source \strong{DEM:}
#' 
#' Ungerechts, L. (2010): DEM 10m (triangulated from aerial photo - b/w). Available online:
#' 
#' \code{http://www.tropicalmountainforest.org/data_pre.do?citid=901}
#' 
#' Jordan, E., Ungerechts, L., Caceres, B. Penafiel, A. and Francou, B. 
#' (2005): Estimation by photogrammetry of the glacier recession on the Cotopaxi Volcano
#' (Ecuador) between 1956 and 1997. \emph{Hydrological Sciences}, 50, 949-961
#' 
#' \strong{Landslide Data:}
#' 
#' Muenchow, J., Brenning, A., Richter, R. (2012) Geomorphic process rates of landslides along a humidity gradient
#' in the tropical Andes, Geomorphology, 139-140, 271-284
#' 
#' Stoyan, R. (2000). Aktivitat, Ursachen und Klassifikation der Rutschungen in San Francisco/Sudecuador.
#' Unpublished Diploma Thesis, University of Erlangen-Nuremberg, Germany.
#' @examples
#' \dontrun{
#' library(RSAGA)
#' data(landslides)
#' 
#' # Print the DEM header:
#' dem$header
#' 
#' # Write the DEM to a SAGA grid:
#' write.sgrd(data = dem, file = "dem", header = dem$header, env = env)
#' 
#' # Calculate slope of DEM:
#' rsaga.slope(in.dem = "dem", out.slope = "slope", method = "poly2zevenbergen", env = env)
#' 
#' # Pick slope values at landslide points,
#' # added to landslides data.frame as variable "slope":
#' landslides <- pick.from.saga.grid(data = landslides,
#'                                   filename = "slope",
#'                                   varname = "slope",
#'                                   env = env)
#' }
NULL