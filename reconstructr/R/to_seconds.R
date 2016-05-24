#'@title to_seconds
#'@description converts POSIX timestamps, or character representations of other timestamp formats, to
#'their numeric value (represented in seconds).
#'
#'@details \code{to_seconds} is designed to enable the rapid conversion of timestamps into their representation
#'as seconds, enabling them to be consumed by \code{\link{reconstruct_sessions}}.
#'
#'@param x a vector of POSIXlt/POSIXct timestamps, or character strings representing a different timestamp
#'format.
#'
#'@param format the format the timestamps take - see \code{\link{strptime}}. Does not need to be set for
#'POSIXlt or POSIXct timestamps.
#'
#'@return a vector of second-values, one for each timestamp.
#'
#'@seealso \code{\link{reconstruct_sessions}} for making use of the new values.
#'
#'@examples
#'#Converting non-POSIX timestamps to seconds
#'data("session_dataset")
#'session_dataset$timestamp <- to_seconds(x = session_dataset$timestamp, format = "%Y%m%d%H%M%S")
#'
#'#Converting POSIX timestamps
#'current_time_in_seconds <- to_seconds(x = Sys.time())
#'@export
to_seconds <- function(x, format){
  if(any(c("POSIXlt","POSIXct") %in% class(x))){
    return(as.numeric(x))
  } else {  
  return(as.numeric(strptime(x,format)))
  }
}