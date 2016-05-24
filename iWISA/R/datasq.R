#' Magnetometer Observations of Station ABG, PHU, TUC and FRD from March to April 2001
#'
#' The H-component of the Ground-based magnetometer records from four stations, Alibag (ABG), Phuthuy (PHU), Tucson (TUC) and  Fredericksburg (FRD), used to estimate Sq variation. 
#' This data frame has 87840 rows, corresponding to a two months interval (March-April 2001), and 4 columns corresponding 
#' to each station.
#'
#'
#' @docType data
#'
#' @usage data(datasq)
#'
#' @format A data frame containing 87840 observations for each four stations.
#'
#' @keywords datasets
#'
#' @source \href{http://www.intermagnet.org/data-donnee/download-eng.php}{INTERMAGNET}
#'
#' @examples
#' data(datasq)
#' head(datasq)
#' # plot of observations in station ABG
#' plot(datasq[,1])
"datasq"