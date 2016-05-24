#' Extract element from a vague date
#' 
#' The NBN Gateway stores dates in three fields: \code{startDate}, 
#' \code{endDate} and \code{dateType}. This allows for storage of date with 
#' various degrees of uncertainty (referred to as "vague dates"). This function
#' allows the \code{year}, \code{month}, \code{week} or \code{day} to be
#' extracted from the vague date whilst properly taking into account the
#' dateType.
#' 
#' @export
#' @param startDate Date returned by NBN Gateway as a string in the format
#'   \code{yyyy-mm-dd}
#' @param endDate Date returned by NBN Gateway as a string in the format 
#'   \code{yyyy-mm-dd}
#' @param dateType One or two letter code, returned by NBN Gateway, describing 
#'   the nature of the vague date (see details)
#' @param part The type of information you want to extract. Possibilities are 
#'   \code{year} (default), \code{month}, \code{week} or \code{day}
#' @return An integer number. -1 is returned if the requested element
#' is not available or does not make sense for the given dateType (see details)  
#' @details
#' The start and end dates are returned by the NBN Gateway as strings in the
#' format \code{yyyy-mm-dd}. The dateType is a one or two letter code as
#' follows:
#' \itemize{
#' \item{\code{D} - \bold{normal date} - start and end dates are the same} 
#' \item{\code{DD} - \bold{date range} (e.g. the period for which a trap was
#' set). startDate is the first day in the range and endDate the last}
#' \item{\code{O} - \bold{month} - startDate is the first of the month, endDate
#' is the last day of the month}
#' \item{\code{OO} - \bold{month range} - startDate is the first day of the
#' starting month of the range, endDate is the last day of the last month in the
#' range}
#' \item{\code{Y} - \bold{year} - startDate is 01 Jan and endDate 31 Dec of the
#' year}
#' \item{\code{YY} - \bold{year range} - startDate is 01 Jan of the year at the
#' start of the range and endDate is 31 Dec of the last year of the range}
#' \item{\code{-Y} - \bold{to year} - represents an uncertain period where only
#' the ending year is known, e.g. the year of publication of a book from which
#' an observation was extracted which gives no indication of when the
#' observation was made - so all we know is that it was made before the book was
#' published. startDate is NULL and endDate is 31 Dec of the year in question.}
#' \item{\code{Y-} - \bold{since year} - represents an uncertain date where all 
#' we know is that it was made after some year. endDate is NULL and startDate is
#' 01 Jan of the year in question. This is \bold{deprecated} (use a year range
#' closed by the year in which the record was extracted (because we know the 
#' observation must have been made before then!), but dates in this format exist
#' on the Gateway.}
#' \item{\code{U} or \code{ND} - \bold{Unknown} - the date was unknown; both
#' start and endDate are NULL}
#' \item{\code{M} - \bold{month} - represents a month when the year is not known. startDate is
#' the first of the month in the year 9999, enddate is the last day of the month
#' in 9999}
#' \item{\code{S} - \bold{season} - represents a seaon when the year is not
#' known. Seasons are taken as a month ranges and can be as follows:
#' \code{spring} - March-May, \code{summer} - June-Aug, \code{autumn} - Sep-Nov
#' or \code{winter} - Dec-Feb. startDate is the first day of the month at the
#' start of the range in year 9999 and endDate is the last day of the month at
#' the end of the range in 9999}
#' \item{\code{P} \bold{publication date} - year of publication - equivalent to
#' \code{Y}}
#' }
#' The value that is returned depends both on the element that is requested and
#' the dateType of the vague date that is passed to the function as follows:
#' \itemize{
#' \item{\code{year} - this is normally the year from endDate. It cannot be
#' determined in the case of \code{Y-}, \code{U}, \code{ND}, \code{S} or
#' \code{M} date types so -1 will be returned in these cases}
#' \item{\code{month} - this is normally the month from endDate. It is only
#' available if the dateType is \code{D}, \code{DD} (and both dates are in the
#' same month), \code{O} or \code{M}. Otherwise -1 is returned.}
#' \item{\code{week} - week number within year. This is normally the week number
#' of the end date. It is only available if the dateType is \code{D} or
#' \code{DD} and both dates are in the same week. Otherwise -1 is returned.}
#' \item{\code{day} - day number within year. This is only available for
#' dateType \code{D}. Otherwise -1 is returned.}
#' }
#' @author Stuart Ball, JNCC \email{stuart.ball@@jncc.gov.uk}
#' @examples
#' datePart("2010-06-15", "2010-06-15", "D", "year") # returns 2010
#' datePart("2010-06-15", "2010-06-15", "D", "mon") # returns 6
#' 
datePart <- function(startDate, endDate, dateType, part=c("year","month","week","day")) {
    
    ##----------------------------------------------------------
    ## get the day number from string date in yyyy-mm-dd format
    dayNumber <- function(date) {
        d <- as.POSIXlt(date)
        d$yday
    }
    
    ## initialise
    part <- match.arg(part)
    v <- -1
    
    ## We can never get anything useful from some date types 
    if (!dateType %in% c("U", "ND", "S")) {
        switch(part,
               year={
                   if (!dateType %in% c("Y-","M")) {
                       v <- as.integer(substr(endDate, 1, 4))
                   }       
               },
               month={
                   if (dateType %in% c("D","DD","O","M")) {
                       s <- as.integer(substr(startDate, 6, 7))
                       e <- as.integer(substr(endDate, 6, 7)) 
                       if (s == e) v <- e
                   }
               },
               week={
                   if (dateType %in% c("D","DD")) {
                       s <- dayNumber(startDate) %/% 7
                       e <- dayNumber(endDate) %/% 7 
                       if (s == e) v <- e
                   }
               },
               day={
                     if (dateType == "D") {
                         v <- dayNumber(endDate)
                     }
               }
               )
    }
    return(v)
}
