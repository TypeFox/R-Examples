#' @include as.precintcon.daily.r
NULL

#' @name as.daily
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @aliases as.precintcon.daily as.daily
#' 
#' @title Converting a data.frame to a daily precipitation serie
#' @description Converts a \code{data.frame} to a \code{precintcon.daily}.
#' @param object a precintcon.daily or a data.frame containing 33 columns.
#' @param na.value the value used for representing non-existent values
#' (Default value: NA).
#' @usage as.daily(object, na.value = NA)
#' @return A \code{data.frame} (\code{precintcon.daily}) containing the 
#' following variables:
#' \itemize{
#' \item \code{year} is the year.
#' \item \code{month} is the month.
#' \item \code{d1} is the precipitation value in millimeters of the 1st day 
#' of the month.
#' \item \code{d2} is the precipitation value in millimeters of the 2nd day 
#' of the month.
#' \item \code{d3} is the precipitation value in millimeters of the 3rd day 
#' of the month.
#' \item \code{d4} is the precipitation value in millimeters of the 4th day 
#' of the month.
#' \item \code{d5} is the precipitation value in millimeters of the 5th day 
#' of the month.
#' \item \code{d6} is the precipitation value in millimeters of the 6th day 
#' of the month.
#' \item \code{d7} is the precipitation value in millimeters of the 7th day 
#' of the month.
#' \item \code{d8} is the precipitation value in millimeters of the 8th day 
#' of the month.
#' \item \code{d9} is the precipitation value in millimeters of the 9th day 
#' of the month.
#' \item \code{d10} is the precipitation value in millimeters of the 10th day 
#' of the month.
#' \item \code{d11} is the precipitation value in millimeters of the 11th day 
#' of the month.
#' \item \code{d12} is the precipitation value in millimeters of the 12th day 
#' of the month.
#' \item \code{d13} is the precipitation value in millimeters of the 13th day 
#' of the month.
#' \item \code{d14} is the precipitation value in millimeters of the 14th day 
#' of the month.
#' \item \code{d15} is the precipitation value in millimeters of the 15th day 
#' of the month.
#' \item \code{d16} is the precipitation value in millimeters of the 16th day 
#' of the month.
#' \item \code{d17} is the precipitation value in millimeters of the 17th day 
#' of the month.
#' \item \code{d18} is the precipitation value in millimeters of the 18th day 
#' of the month.
#' \item \code{d19} is the precipitation value in millimeters of the 19th day 
#' of the month.
#' \item \code{d20} is the precipitation value in millimeters of the 20th day 
#' of the month.
#' \item \code{d21} is the precipitation value in millimeters of the 21th day 
#' of the month.
#' \item \code{d22} is the precipitation value in millimeters of the 22th day 
#' of the month.
#' \item \code{d23} is the precipitation value in millimeters of the 23th day 
#' of the month.
#' \item \code{d24} is the precipitation value in millimeters of the 24th day 
#' of the month.
#' \item \code{d25} is the precipitation value in millimeters of the 25th day 
#' of the month.
#' \item \code{d26} is the precipitation value in millimeters of the 26th day 
#' of the month
#' \item \code{d27} is the precipitation value in millimeters of the 27th day 
#' of the month.
#' \item \code{d28} is the precipitation value in millimeters of the 28th day 
#' of the month.
#' \item \code{d29} is the precipitation value in millimeters of the 29th day 
#' of the month.
#' \item \code{d30} is the precipitation value in millimeters of the 30th day 
#' of the month.
#' \item \code{d31} is the precipitation value in millimeters of the 31th day 
#' of the month. 
#' }
#' @seealso 
#' \code{\link{as.decade}}
#' \code{\link{as.annual}}
#' \code{\link{as.seasonal}}
#' \code{\link{as.monthly}} 
#' @usage as.daily(object, na.value = NA)
#' @examples
#' ##
#' # Loading the daily precipitation serie.
#' data(daily)
#' 
#' ##
#' # Converting precipitation
#' as.daily(daily)
#' @keywords daily precipitation
#' @export
as.daily <- function(object, na.value = NA) {
   return(as.precintcon.daily(object, na.value))
}