#' Atmospheric data.
#'
#' @format
#' \describe{
#' \item{lat,long}{Location of measurement. Evenly spaced 24 by 24 spatial
#'   grid from longitude 113.8W to 56.2W and from latitude 36.2N to 21.2S}
#' \item{year,month}{72 points in time; once per month from Jan 1995 to
#'   Dec 2000.}
#' \item{surftemp}{Mean Surface Temperature from Clear Sky Composite (ts):
#'   The monthly mean temperature based on the energy being emitted from the
#'   Earth's surface under clear sky conditions (in K).}
#' \item{temp}{Mean Near-Surface Air Temperature (tsa_tovs):
#'   The monthly mean temperature of the air near the surface of the Earth
#'   (in K).}
#' \item{pressure}{Mean Surface Pressure (ps_tovs): The monthly mean
#'   atmospheric surface pressure at a given location on the Earth's surface.
#'   (in mb)}
#' \item{ozone}{Mean Ozone Abundance (o3_tovs): The monthly mean amount of
#'   total ozone in the atmospheric column (in Dobsons)}
#' \item{cloudlow}{Mean Low Cloud Coverage (ca_low): The monthly mean percent
#'   of the sky covered by clouds with cloud top pressure greater than 680 mb
#'   or roughly less than 3.24 km.}
#' \item{cloudmid}{Mean Medium Cloud Coverage (ca_med): The monthly mean
#'   percent of the sky covered by clouds with cloud top pressure between 440 -
#'   680 mb or roughly 3.24 to 6.5 km.}
#' \item{cloudhigh}{Mean High Cloud Coverage (ca_high): The monthly mean
#'   percent of the sky covered by clouds with cloud top pressure less than or
#'   equal to 440 mb or roughly greater than 6.5 km.}
#' }
#' @source \url{http://stat-computing.org/dataexpo/2006/}
"atmos"

#' Elevation.
#'
#' As of Jan 1998.
#'
#' @format
#' \describe{
#' \item{lat,long}{Location of measurement. Evenly spaced 24 by 24 spatial
#'   grid from longitude 113.8W to 56.2W and from latitude 36.2N to 21.2S}
#' \item{elev}{Height of location above reference plane (in m)}
#' }
#' @source \url{http://stat-computing.org/dataexpo/2006/}
"elev"

