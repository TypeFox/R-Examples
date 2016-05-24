#' A Timing Function for SPMD Routines
#' 
#' A timing function for use with parallel codes executed in the batch SPMD
#' style.
#' 
#' Finds the min, mean, and max execution time across all independent processes
#' executing the operation \code{timed}.
#' 
#' @param timed
#' expression to be timed.
#' 
#' @return
#' A named vector containing the minimum, mean, and maximum time across
#' all processors in the communicator.  All values are global.
#' 
#' @keywords Timing
#' @export
timer <- function(timed)
{
  ltime <- system.time(timed)[3]
  
  mintime <- allreduce(ltime, op='min')
  maxtime <- allreduce(ltime, op='max')
  
  meantime <- allreduce(ltime, op='sum') / comm.size()
  
  return( c(min=mintime, mean=meantime, max=maxtime) )
}
