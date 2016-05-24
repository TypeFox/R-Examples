
#' @title Retrieve Data from the HYDRA Project
#'
#' @description Retrieves rainfall data from the City of Portland's HYDRA monitoring stations.
#'
#' @seealso \code{\link{hail_hydra}} for data retrieval.
#'
#' @docType package
#' @name hail
NULL

#' @title HYDRA metadata
#'
#' @description metadata about the HYDRA rain gages.
#'
#' @format A data frame with 53 rows and 3 variables:
#' \describe{
#'   \item{station}{The station name. This can be passed into \code{\link{hail_hydra}} to
#'   retrieve data from that station.}
#'   \item{address}{The station's address/}
#'   \item{url}{The URL of the station's raw, unread data/}
#' }
#' @source \url{http://or.water.usgs.gov/precip/}
"hydra_data"