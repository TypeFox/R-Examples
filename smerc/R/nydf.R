#' @name nydf
#' @title Leukemia data for 281 regions in New York.
#' @description This data set contains 281 observations related to leukeumia cases in an 8 county area of the state of New York.  The data were made available in Waller and Gotway (2005) and details are provided there.  These data are related to a similar data set in Waller et al. (1994).  The longitude and latitude coordinates are taken from the NYleukemia data set in the SpatialEpi package for plotting purposes.
#' 
#' @docType data
#' @usage data(nydf)
#' 
#' @format A data frame with 281 rows and 4 columns:
#' \describe{
#'  \item{longitude}{The longitude of the region centroid.  These are NOT the original values provided by Waller and Gotway (2005), but are the right ones for plotting correctly.}
#'  \item{latitude}{The latitude of the region centroid.  These are NOT the original values provided by Waller and Gotway (2005), but are the right ones for plotting correctly.}
#'  \item{population}{The population (1980 census) of the region.}
#'  \item{cases}{The number of leukemia cases between 1978-1982.}
#'  \item{x}{The original 'longitude' coordinate provided by Waller and Gotway (2005).}
#'  \item{y}{The original 'latitude coordinate provided by Waller and Gotway (2005).}
#' }
#' @source Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial Statistics for Public Health Data.  Hoboken, NJ: Wiley.
#' @references Waller, L.A., Turnbull, B.W., Clark, L.C., and Nasca, P. (1994) "Spatial Pattern Analysis to Detect Rare Disease Clusters" in Case Studies in Biometry, N. Lange, L. Ryan, L. Billard, D. Brillinger, L. Conquest, and J. Greenhouse (eds.) New York: John Wiley and Sons.
NULL