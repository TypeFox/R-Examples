#' Houston flights data
#'
#' This dataset contains all flights departing from Houston airports IAH
#' (George Bush Intercontinental) and HOU (Houston Hobby). The data comes
#' from the Research and Innovation Technology Administration at the
#' Bureau of Transporation statistics:
#' \url{http://www.transtats.bts.gov/DatabaseInfo.asp?DB_ID=120&Link=0}
#'
#' \code{src_hflights} caches a SQLite version of the data in a standard
#' location for use in examples.
#'
#' @section Variables:
#'
#' \itemize{
#'   \item \code{Year}, \code{Month}, \code{DayofMonth}: date of departure
#'   \item \code{DayOfWeek}: day of week of departure (useful for removing
#'     weekend effects)
#'  \item \code{DepTime}, \code{ArrTime}: departure and arrival times
#'    (in local time, hhmm)
#'  \item \code{UniqueCarrier}: unique abbreviation for a carrier
#'  \item \code{FlightNum}: flight number
#'  \item \code{TailNum}: airplane tail number
#'  \item \code{ActualElapsedTime}: elapsed time of flight, in minutes
#'  \item \code{AirTime}: flight time, in minutes
#'  \item \code{ArrDelay}, \code{DepDelay}: arrival and departure delays,
#'    in minutes
#'  \item \code{Origin}, \code{Dest} origin and destination airport codes
#'  \item \code{Distance}: distance of flight, in miles
#'  \item \code{TaxiIn}, \code{TaxiOut}: taxi in and out times in minutes
#'  \item \code{Cancelled}: cancelled indicator: 1 = Yes, 0 = No
#'  \item \code{CancellationCode}: reason for cancellation: A = carrier,
#'    B = weather, C = national air system, D = security
#'  \item \code{Diverted}: diverted indicator: 1 = Yes, 0 = No
#' }
#' @docType data
#' @name hflights
#' @usage hflights
#' @format A data frame with 227,496 rows and 21 columns.
#' @examples
#' head(hflights)
NULL

