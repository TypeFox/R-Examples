#' This function retrieves Catchment Mean Rainfall (CMR).
#'
#' @author Claudia Vitolo
#'
#' @description Given the station ID number(s), this function retrieves data (time series in zoo format with accompanying metadata) from the WaterML2 service on the NRFA database. Catchment Mean Rainfall is measured in mm/month.
#'
#' @param id station ID number(s), each number should be in the range [3002,236051].
#' @param metadata Logical, FALSE by default. If metadata = TRUE means that the result for a single station is a list with two elements: data (the time series) and meta (metadata).
#' @param parallel Logical, FALSE by default. If parallel = TRUE means that the function can be used in parallel computations.
#'
#' @return list composed of as many objects as in the list of station ID numbers. Each object can be accessed using their names or index (e.g. x[[1]], x[[2]], and so forth). Each object contains a zoo time series.
#'
#' @examples
#' CMR(18019)
#' # CMR(c(54022,54090,54091))
#'

CMR <- function(id, metadata = FALSE, parallel = FALSE){

  # require(RCurl)
  # require(XML2R)
  # require(stringr)
  # require(zoo)
  # id <- c(54022,54090,54091)

  rain <- getTS(id, type = "cmr", metadata, parallel)

  return(rain)

}
