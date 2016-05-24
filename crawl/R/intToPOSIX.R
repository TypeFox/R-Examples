#' Reverse as.numeric command that is performed on a vector of type POSIXct
#' 
#' Takes integer value produced by \code{as.numeric(x)}, where \code{x} is a
#' POSIXct vector and returns it to a POSIXct vector
#' 
#' @aliases intToPOSIX
#' @param timeVector A vector of integers produced by as.numeric applied to a
#' PSIXct vector
#' @param tz Time zone of the vector (see \code{\link{as.POSIXct}}).
#' @return POSIXct vector
#' @note There is no check that as.numeric applied to a POSIX vector produced
#' \code{timeVector}.  So, caution is required in using this function. It was
#' included simply because I have found it useful
#' @author Devin S. Johnson
#' @examples
#' 
#' #library(crawl)
#' timeVector <- as.numeric(Sys.time())
#' timeVector
#' intToPOSIX(timeVector, tz="")
#' @export
`intToPOSIX` <-
function(timeVector, tz='GMT')
################################################################################
# Convert integer time to POSIX
# timeVector  = integer time in seconds
# tz          = time zone code
################################################################################
{
  Epoch = as.POSIXct(strptime("1970-01-01 00:00:00", "%Y-%m-%d %H:%M:%S",tz=tz),tz=tz)
  Epoch+timeVector
}

