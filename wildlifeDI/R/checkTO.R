# ---- roxygen documentation ----
#
#' @title Check for temporal overlap
#'
#' @description
#' The function \code{checkTO} is a simple function for identifying if, and for how long, two telemetry datasets overlap (temporally) with each other. The function returns a list with three pieces of information: a logical variable indicating if the two trajectories overlap temporally, and timings of the beginning and end of the overlap period.  
#' 
#' @details
#' The function \code{checkTO} can be used to identify if, when, and for how long two telemetry datasets overlap temporally.   
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#'
#' @return
#' A \code{list} of with three pieces of information, whether the two trajectories overlap (\code{$TO}) a logical vector, the beginning (\code{$TOstart}), and end (\code{$TOend}) of the overlap period, stored as \code{POSIX} objects.  
#'
#' @keywords temporal overlap
#' @seealso GetSimultaneous
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' spts <- checkTO(deer37, deer38)
#' 
#' @export
#

checkTO <- function(traj1,traj2){
  #convert to dataframe
  tr1 <- ld(traj1)
  tr2 <- ld(traj2)
  
  #identify the temporal overlap points.
  t.min <- max(c(min(tr1$date),min(tr2$date)))
  t.max <- min(c(max(tr1$date),max(tr2$date)))
  #check to see if there is an overlap period
  if (t.min > t.max){
    TO <- FALSE
    t.min <- t.max <- NA
  } else {
    TO <- TRUE
  }
  
  return(list(TO=TO, TOstart=t.min, TOend=t.max))
}