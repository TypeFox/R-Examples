#' Dow Jones Index Data
#'
#' @format  A data frame with 252 observations on the following 7 variables containing data
#' from  2014-01-01 to 2015-01-01. 
#'  \describe{
#'    \item{\code{Date}}{Date of observation in character string format "\%m/\%d/\%Y" }
#'    \item{\code{DJI.Open}}{Opening value of DJI on the specified date}
#'    \item{\code{DJI.High}}{High value of the DJI on the specified date}
#'    \item{\code{DJI.Low}}{Low value of the DJI on the specified date}
#'    \item{\code{DJI.Close}}{Closing value of the DJI on the specified date}
#'    \item{\code{DJI.Volume}}{the number of shares or contracts traded}
#'    \item{\code{DJI.Adj.Close}}{Close price adjusted for dividends and splits}
#'  }
#' @source Data obtained using \code{yahooSeries} from the \code{fImport} package.
"dowJonesData"