# Data documentation

#ROxygen comments ----
#' Simulated nonstationary dataset
#'
#' A data set containing the necessary components to fit the nonstationary
#' spatial model, simulated from the true model.
#'
#' @format A list with the following objects:
#' \describe{
#'    \item{sim.locations}{A matrix of longitude/latitude coordinates of the
#'    simulated locations.}
#'    \item{mc.locations}{A matrix of longitude/latitude coordinates of the
#'    mixture component locations.}
#'    \item{mc.kernel}{A three-dimensional array, containing the true 2 x 2
#'    kernel covariance matrices for each mixture component location.}
#'    \item{kernel.ellipses}{A three-dimensional array, containing the true 2 x 2
#'    kernel covariance matrices for each simulated location.}
#'    \item{sim.data}{A matrix of the simulated data; each of the ten columns
#'    correspond to an independent and identically distribured replicate.}
#'    \item{lambda.w}{Scalar; the value of the tuning parameter used in the
#'    weight function.}
#'    \item{holdout.index}{Vector; indicates which of the simulated locations
#'    should be used in the hold-out sample.}
#'
#' }
"simdata"


#' Annual precipitation measurements from the western United States, 1997
#'
#' A data set containing the annual precipitation for 1270 locations
#' in the western United States.
#'
#' @format A data frame with the following variables:
#' \describe{
#'    \item{longitude}{Longitude of the monitoring site.}
#'    \item{latitude}{Latitude of the monitoring site.}
#'    \item{annual.ppt}{Annual precipitation for the monitoring site, in
#'    millimeters.}
#'    \item{log.annual.ppt}{Annual precipitation for the monitoring site, in
#'    log millimeters.}
#'
#' }
#' @source \url{http://www.image.ucar.edu/GSP/Data/US.monthly.met/}
"USprecip97"

#' Prediction locations for the western United States
#'
#' A matrix with two columns containing a fine grid of locations for
#' which to make a filled-in prediction map for the western United States.
#'
#' @format A matrix with two columns:
#' \describe{
#'    \item{Column 1}{Longitude of the prediction grid.}
#'    \item{Column 2}{Latitude of the prediction grid.}
#'
#' }
"US.prediction.locs"

#' Mixture component grids for the western United States
#'
#' A list of two mixture component grids for fitting the nonstationary
#' model to the western United States precipitation data.
#'
#' @format A list with two elements:
#' \describe{
#'    \item{Element 1}{Coarse mixture component grid.}
#'    \item{Element 2}{Fine mixture component grid.}
#'
#' }
"US.mc.grids"






