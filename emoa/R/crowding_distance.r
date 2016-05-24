##
## crowding_distance.r - calculate crowding distance.
##

##' Calculate crowding distances.
##'
##' @aliases crowding_distance
##' @title Crowding Distance
##' @export 
##' @param front matrix of function values.
##'
##' @return crowding distance for each function value.
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
crowding_distance <- function(front)
  .Call(do_crowding_distance, front)
