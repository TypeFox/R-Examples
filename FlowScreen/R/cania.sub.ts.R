#' @name cania.sub.ts
#' @title Subset of the Caniapiscau River Daily Flows
#' @description This data set includes a subset of the mean daily streamflow 
#' for the Caniapiscau Rivers. It includes observations from 1970-1995
#' (hydrologic years). The code used to subset and modify the original
#' data is shown below. 
#' @docType data
#' @usage data(caniapiscau)
#' @format Formatted as a data.frame with the following columns:
#'   \itemize{
#'     \item ID - Water Survey Canada Station ID
#'     \item Date - Date of observation, formatted as YYYY-mm-dd
#'     \item Flow - Mean daily streamflow, measured in m3/s
#'     \item Code - Data Quality Code
#'     \item Agency - Source Agency (Water Survey Canada)
#'     \item Year - Calendar year
#'     \item month - Calendar month
#'     \item doy - Calendar day of year
#'     \item hyear - Hydrologic year
#'     \item hmonth - Hydrologic month
#'     \item hdoy - Hydrologic day of year
#'   }
#' @source Environment Canada. 2010. EC Data Explorer V1.2.30. \cr
#'   Water Survey of Canada V1.2.30 https://www.ec.gc.ca/rhc-wsc/
#' @examples
#' # Code used to subset and modify original Caniapiscau series:
#' \dontrun{
#' data(caniapiscau)
#' cania.ts <- create.ts(caniapiscau, hyrstart=3)
#' cania.sub.ts <- subset(cania.ts, cania.ts$hyear %in% c(1970:1995))
#' }
#' # example use of example subset flow series
#' data(cania.sub.ts)
#' head(cania.sub.ts)
#' str(cania.sub.ts)

"cania.sub.ts"