#' @name caniapiscau
#' @title Caniapiscau River Daily Flows
#' @description This data set includes the mean daily streamflow for the Caniapiscau
#'   River. The file has been read from it's original .csv format using  
#'   \code{\link{read.flows}}. The Caniapiscau River is located in Nunavik, Quebec,
#'   Canada, and flows northward. The headwaters (representing 45 percent of the 
#'   total flow) were dammed to create the Caniapiscau Reservoir, which started 
#'   filling in 1981. In 1985, the reservoir was diverted to the west into the 
#'   La Grande hydroelectric complex.
#'   This flow time series is used as an example of a river with a known 
#'   change point to demonstrate the package's screening capabilities.
#' @docType data
#' @usage data(caniapiscau)
#' @format Formatted as a data.frame with the following columns:
#'   \itemize{
#'     \item ID - Water Survey Canada Station ID
#'     \item PARAM - Parameter ID (1 indicates flow)
#'     \item Date - Date of observation, formatted as YYYY-mm-dd
#'     \item Flow - Mean daily streamflow, measured in m3/s
#'     \item Agency - Source Agency (Water Survey Canada)
#'   }
#' @source Environment Canada. 2010. EC Data Explorer V1.2.30. \cr
#'   Water Survey of Canada V1.2.30 https://www.ec.gc.ca/rhc-wsc/
#' @examples
#' data(caniapiscau)
#' head(caniapiscau)
#' str(caniapiscau)

"caniapiscau"