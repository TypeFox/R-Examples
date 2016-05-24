FitHTTPHeader <- function(token){

  return(c(Authorization = paste("Bearer",token,sep = ' ')))

  }

ValidateDateRange <- function(startDate, endDate) {

  stopifnot(inherits(startDate,"POSIXct"))

  stopifnot(inherits(endDate,"POSIXct"))

}

ValidateDatasource <- function(datasource) {
  
  stopifnot(length(datasource)==1)
  
}

EpochTime <- function(inDate) {
  
  return (as.integer64(as.numeric(inDate) * 1000000000))
  #return (paste(as.numeric(inDate),"000000000",sep = ""))
  
}

#' @title NanosToPOSIXct
#' @rdname NanosToPOSIXct
#' @export
#'
#' @param nanos - Nanoseconds from epoch
#' @description
#' Converts nanoseconds from epoch (as provided by Google Fit) to POSIXct
#' @examples
#' NanosToPOSIXct(1388534400000000000)

NanosToPOSIXct <- function(nanos) {
  
  return (as.POSIXct(round(as.integer64(nanos) /
                              as.integer64(1000000000)),
                     origin = "1970-01-01")
  )
  
}
  