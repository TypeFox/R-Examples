#' @name co
#' @title Geochemical measurements for 960 sites in Colorado.
#' @description These data were collected by the United States Geological Survey (USGS).  Their description is as follows:  In 2006, soil samples were collected at 960 sites (1 site per 280 square kilometers) throughout the state of Colorado. These samples were collected from a depth of 0 to 15 centimeters and, following a near-total multi-acid digestion, were analyzed for a suite of more than 40 major and trace elements. The resulting data set provides a baseline for the natural variation in soil geochemistry for Colorado and forms the basis for detecting changes in soil composition that might result from natural processes or anthropogenic activities. 
#' 
#' Latitude and Longitude determined by hand-held GPS instrument using WGS84 datum.  The longitude and latitude coordinates were converted to UTM coordinates using the following commands from the \code{sp} package:
#' lonlat = co[, c("longitude", "latitude")]; 
#' coordinates(lonlat) <- c("longitude", "latitude")
#' proj4string(lonlat) <- CRS("+proj=longlat +datum=WGS84")  ## for example
#' xy <- spTransform(lonlat, CRS("+proj=utm +zone=13 ellps=WGS84"))
#' 
#' @docType data
#' @usage data(co)
#' 
#' @format A data frame with 960 rows and 31 columns:
#' \describe{
#'  \item{easting}{Easting (m)}
#'  \item{northing}{Northing (m)}
#'  \item{latitude}{The latitude of the site.}
#'  \item{longitude}{The longitude of the site.}
#'  \item{Al}{aluminum (percent)}
#'  \item{Ca}{calcium (percent)}
#'  \item{Fe}{iron (percent)}
#'  \item{K}{potassium (percent)}
#'  \item{Mg}{magnesium (percent)}
#'  \item{Na}{sodium (percent)}
#'  \item{Ti}{titanium (percent)}
#'  \item{Be}{beryllium (mg/kg)}
#'  \item{Ce}{cerium (mg/kg)}
#'  \item{Co}{cobalt (mg/kg)}
#'  \item{Cr}{chromium (mg/kg)}
#'  \item{Cu}{copper (mg/kg)}
#'  \item{Ga}{gallium (mg/kg)}
#'  \item{La}{lanthanum (mg/kg)}
#'  \item{Li}{lithium (mg/kg)}
#'  \item{Mo}{manganese (mg/kg)}
#'  \item{Nb}{molybdenum (mg/kg)}
#'  \item{Ni}{niobium (mg/kg)}
#'  \item{Rb}{rubidium (mg/kg)}
#'  \item{Sc}{scandium (mg/kg)}
#'  \item{Sn}{tin (mg/kg)}
#'  \item{Th}{thorium (mg/kg)}
#'  \item{Tl}{thallium (mg/kg)}
#'  \item{U}{uranium (mg/kg)}
#'  \item{V}{vanadium (mg/kg)}
#'  \item{W}{tungsten (mg/kg)}
#'  \item{Y}{yttrium (mg/kg)}
#' }
#' @source U.S. Geological Survey, Data Series 520, 9 p.  \url{http://pubs.usgs.gov/ds/520/}.
#' @references Smith, D.B., Ellefsen, K.J., and Kilburn, J.E., 2010, Geochemical data for Colorado soils -- Results from the 2006 state-scale geochemical survey: U.S. Geological Survey, Data Series 520, 9p.
NULL