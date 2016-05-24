#' Magnetometer Observations of Station HER, KAK, HON and SJG from March to April 2001
#'
#' The H-component of the Ground-based magnetometer records four stations,Hermanus (HER), Kakioka (KAK), Honolulu (HON) and San Juan (SJG), used to 
#' estimate WISA and preindex. This data frame has 87840 rows, corresponding to a 
#' two months (March to April 2001) of 1-minute data, and 4 columns corresponding to each station.
#'
#' @docType data
#'
#' @usage data(record)
#'
#' @format A data frame containing 87840 observations of H-component for four stations.
#'
#' @keywords datasets
#'
#' @source \href{http://www.intermagnet.org/data-donnee/download-eng.php}{INTERMAGNET}
#' 
#' @examples
#' data(record)
#' head(record)
#' 
#' #plot of observations from station HER
#' plot(record[,1])
"record"