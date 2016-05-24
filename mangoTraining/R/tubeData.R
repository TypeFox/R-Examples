#' London Tube Performace data
#' 
#' @format  A data frame with 1050 observations on the following 9 variables.
#'  \describe{
#'    \item{\code{Line}}{A factor with 10 levels, one for each London tube line}
#'    \item{\code{Month}}{A numeric vector indicating the month of the observation}
#'    \item{\code{Scheduled}}{A numeric vector giving the scheduled running time}
#'    \item{\code{Excess}}{A numeric vector giving the excess running time}
#'    \item{\code{TOTAL}}{A numeric vector giving the total running time}
#'    \item{\code{Opened}}{A numeric vector giving the year the line opened}
#'    \item{\code{Length}}{A numeric vector giving the line length}
#'    \item{\code{Type}}{A factor indicating the type of tube line}
#'    \item{\code{Stations}}{A numeric vector giving the number of stations on the line}
#'  }
#'  
#' @source This data was taken from "http://data.london.gov.uk/datafiles/transport/assembly-tube-performance.xls"
"tubeData"