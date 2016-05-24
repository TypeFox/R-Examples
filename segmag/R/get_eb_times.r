#' Retrieve event boundary times from a segmag object
#' 
#' Returns the times of event boundaries from a segmag object. Event boundaries
#' are defined as the center of local maxima in segmentation magnitude that are
#' higher than a critical cutoff value. The critical cutoff value should be
#' determined with the \code{\link{bootstrap_critical_cutoffs}} function.
#' 
#' @param segmag object of class \code{\link{segmag}}
#' @param cutoff numeric value determining the critical cutoff in segmentation
#'               magnitude
#' @return numeric vector with event boundary times 
#' @seealso \code{\link{get_eb_times_segmag_diff}}, \code{\link{bootstrap_critical_cutoffs}} 
#' @examples #see ?segmag for an example
#' @export
get_eb_times <- function(segmag, cutoff)
{
  if (! is.segmag(segmag)) stop("segmag must be an object of class segmag")
  
  if (! is.numeric(cutoff) || length(cutoff) != 1)
    stop("cutoff must be numeric of length 1")
    
  segmag_df <- segmag$data[flag_maxima_positions(segmag$data$segmentation_magnitude),]
  
  return( segmag_df$time[segmag_df$segmentation_magnitude >= cutoff] )
}

#' Retrieve event boundary times for a difference of segmag objects
#' 
#' Specific function when calculating event boundaries for the difference of
#' two segmag objects (segmag - segmag_substract). Returns the times of event
#' boundaries that are defined as the center of local maxima (minima) in
#' segmentation magnitude that are higher (lower) than a critical cutoff max/min values.
#' The critical cutoff max/min values should be determined with the
#' \code{\link{bootstrap_critical_cutoffs}} function.
#' 
#' @param segmag object of class \code{\link{segmag}}
#' @param segmag_substract object of class \code{\link{segmag}}.
#' @param cutoff_max numeric value determining the critical cutoff for maxima in
#'                   segmentation magnitude
#' @param cutoff_min numeric value determining the critical cutoff for minima in
#'                   segmentation magnitude
#' @return numeric vector with event boundary times. If cutoff_max or cutoff_min is
#'         NULL than the respective event boundaries are omitted.
#' @seealso \code{\link{get_eb_times}}, \code{\link{bootstrap_critical_cutoffs}} 
#' @export
get_eb_times_segmag_diff <- function(segmag, segmag_substract, cutoff_max = NULL, cutoff_min = NULL)
{
  if (! is.segmag(segmag)) stop("segmag must be an object of class segmag")
  if (! is.segmag(segmag_substract)) stop("segmag_substract must be an object of class segmag")
  
  if (! is.null(cutoff_max) && (! is.numeric(cutoff_max) || length(cutoff_max) != 1))
    stop("cutoff_max must be numeric of length 1 or NULL")
  
  if (! is.null(cutoff_min) && (! is.numeric(cutoff_min) || length(cutoff_min) != 1))
    stop("cutoff_min must be numeric of length 1 or NULL")
  
  if (segmag$time_min != segmag_substract$time_min || segmag$time_max != segmag_substract$time_max || segmag$time_steps != segmag_substract$time_steps) stop("time_min, time_max, and time_steps of segmag and segmag_substract must be equal")
  if (length(segmag$data$time) != length(segmag_substract$data$time) || sum(segmag$data$time != segmag_substract$data$time) > 0) stop("segmag$data$time and segmag_substract$data$time must be equal")
  segmag$data$segmentation_magnitude <- segmag$data$segmentation_magnitude - segmag_substract$data$segmentation_magnitude
  
  segmag_df_maxima <- segmag$data[flag_maxima_positions(segmag$data$segmentation_magnitude),]
  segmag_df_minima <- segmag$data[flag_minima_positions(segmag$data$segmentation_magnitude),]
  
  eb_times <- numeric(0)
  if (! is.null(cutoff_max)) eb_times <- c(eb_times, segmag_df_maxima$time[segmag_df_maxima$segmentation_magnitude >= cutoff_max])
  if (! is.null(cutoff_min)) eb_times <- c(eb_times, segmag_df_minima$time[segmag_df_minima$segmentation_magnitude <= cutoff_min])
  
  return( eb_times )
}