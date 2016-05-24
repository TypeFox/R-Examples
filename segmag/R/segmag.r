#' @useDynLib segmag
#' @importFrom Rcpp evalCpp
NULL

#' Create Segmentation Object
#' 
#' This function creates a segmag object from a vector of participant ids and
#' a vector of times when the participants pressed the segmentation key on the keyboard.
#' All functions in the segmag package work on this object (e.g., plotting results, 
#' determining event boundaries). The additional parameters define the size and
#' offset of the Gaussian that is centered around each key press and the time window
#' in the data set to consider (time_min, time_max, time_steps).
#' 
#' First, segmentation magnitude for each participant across time is calculated 
#' by centering a Gaussian around each key press. If multiple Gaussians overlap
#' across time then only the maximum values is used (not sum) to ensure an equal
#' weight of each participant on the overall segmentation magnitude. Thereafter,
#' the segmentation magnitudes of the participants are summed up to define the
#' overall segmentation magnitude across time. The higher the segmentation magnitude
#' at one point in time the more participants pressed a key around this time
#' point. To account for the fact that participants have a certain temporal error
#' in their key presses, Gaussians are used to expand the influence of a single key press
#' into time. Furthermore, an offset to these Gaussians can be defined in order to
#' account for manual reaction times and to get a better estimate of the "real"
#' time point of an event boundary.
#' 
#' In order to achieve a decent calculation speed, a fixed time scale with interval length
#' time_steps and starting from time_min is used. All key-press times are rounded to their closest
#' interval. A warning is issued if this changes the raw key-press times.
#' 
#' @param ids factor assigning a participant id to each key press (same length as
#'            time_keypresses)
#' @param time_keypresses numeric vector with the times when the segmentation key
#'                        was pressed on the keyboard
#' @param data optional data.frame where ids and time_keypresses are stored in
#' @param time_min time window of key press times used in calculations
#' @param time_max time window of key press times used in calculations
#' @param time_steps interval length of time steps within time window
#' @param gauss_offset offset the overlaid Gaussian relative to the key press times
#'                     in order to account for manual reaction times (0: centered,
#'                     negative values: assume that event boundary occurred before key press)
#' @param gauss_sd sd of overlaid Gaussian
#' @param gauss_cutoff speed up calculations by not considering values of the overlaid Gaussian
#'                     that are more than gauss_cutoff units from the center of the Gaussian away
#'                     (because they are very close to 0)
#' @return segmag object (also contains a data.frame with segmentation magnitude across time as $data)
#' @seealso \code{\link{bootstrap_critical_cutoffs}}, \code{\link{get_eb_times}}, \code{\link{plot.segmag}}
#' @examples
#' # segmentation responses (key presses) of 6 participants watching a movie (30 seconds long)
#' participant_ids <- factor(c(1,1,1,1,2,2,3,3,3,3,4,4,4,5,5,6,6,6))
#' time_keypresses <- c(7,12,18,25,12.1,24.9,6.9,10,25.2,29,7.2,12.05,17.5,7.05,25,6.9,16.1,25)
#' 
#' # create segmag object
#' segmag1 <- segmag(participant_ids, time_keypresses, time_min=0, time_max=30)
#' 
#' \dontrun{
#' # estimate the critical cutoff against an alpha level of 0.05
#' # Note: This is an estimate an will vary slightly against multiple calls of this function
#' #       (variation is the lower the higher n_bootstrap is set)
#' critical_cutoff <- bootstrap_critical_cutoffs(segmag1, 5000, .95)
#' }
#' \dontshow{
#' critical_cutoff <- 1.623081
#' }
#' 
#' # timestamps of significant event boundaries within the movie
#' eb_times <- get_eb_times(segmag1, critical_cutoff)
#' 
#' plot(segmag1, critical_cutoff, eb_times)
#' 
#' @export
segmag <- function(ids, time_keypresses, data = NULL, time_min = (min(min(time_keypresses),(min(time_keypresses)+gauss_offset))-gauss_cutoff), time_max = (max(max(time_keypresses),(max(time_keypresses)+gauss_offset))+gauss_cutoff), time_steps = 0.01, gauss_offset = 0, gauss_sd = 1, gauss_cutoff = 6 * gauss_sd)
{
  # Überlegen.. wie sicherstellen dass bei bootstrap nicht keypresses an zeiten platziert werden, die dann verloren gehen, da über cutoff verschoben werden?
  # indexes calculated in c notation! 0..i-1
  #gauss_cutoff: # computational cutoff because almost 0.. einseitig, also je cutoff nach links und cutoff nach rechts
  
  if (! is.null(data))
  {
    if (! is.data.frame(data)) stop("data must be a data.frame or NULL")
    ids = eval(substitute(ids), data, parent.frame())
    time_keypresses = eval(substitute(time_keypresses), data, parent.frame())
  }
  
  if (is.null(ids) || ! is.factor(ids) || length(ids) == 0) stop("ids must be a factor with length > 0")
  if (is.null(time_keypresses) || ! is.numeric(time_keypresses) || length(time_keypresses) == 0) stop("time_keypresses must be a numeric with length > 0")
  if (length(ids) != length(time_keypresses)) stop("ids and time_keypresses must have the same length")
  
  if (! is.numeric(gauss_offset) || length(gauss_offset) != 1) stop("gauss_offset must be numeric and of length 1")
  if (! is.numeric(gauss_sd) || length(gauss_sd) != 1) stop("gauss_sd must be numeric and of length 1")
  if (! is.numeric(gauss_cutoff) || length(gauss_cutoff) != 1) stop("gauss_cutoff must be numeric and of length 1")
  
  if (! gauss_sd > 0) stop("gauss_sd must be greater than 0")
  if (! gauss_cutoff > 0) stop("gauss_cutoff must be greater than 0")
  
  if (! is.numeric(time_min) || length(time_min) != 1) stop("time_min must be numeric and of length 1")
  if (! is.numeric(time_max) || length(time_max) != 1) stop("time_max must be numeric and of length 1")
  if (! is.numeric(time_steps) || length(time_steps) != 1) stop("time_steps must be numeric and of length 1")
  
  if (! (time_max > time_min)) stop("time_max must be greater than time_min")
  if (! time_steps > 0) stop("time_steps must be greater than 0")
  
  
  index_time_max <- round((time_max - time_min) / time_steps)
  
  if (index_time_max != ( (time_max - time_min) / time_steps ) )
  {
    time_max = time_min + (index_time_max * time_steps)
    warning(paste0("time_max was changed to ",time_max," such that it is a multiple of time_min + n * time_steps"))
  }
  
  indexes_gauss_offset <- round(gauss_offset / time_steps)
  
  if (indexes_gauss_offset != ( gauss_offset / time_steps ) )
  {
    warning(paste0("gauss_offset was changed to ",(indexes_gauss_offset * time_steps)," such that it is a multiple of time_steps"))
  }
  
  # Gauss Array geht von 0 (linker Rand) über 0 + gauss_n_indexes_per_side (Maximum) bis 0 + gauss_n_indexes_per_side + gauss_n_indexes_per_side (rechter Rand)
  # Gauss Array ist also gauss_n_indexes_per_side * 2 + 1 lang
  gauss_n_indexes_per_side <- ceiling(gauss_cutoff / time_steps)
  
  if (gauss_n_indexes_per_side != ( gauss_cutoff / time_steps ) )
  {
    warning(paste0("gauss_cutoff was rounded up to ",(gauss_n_indexes_per_side * time_steps)," such that it is a multiple of time_steps"))
  }
  
  gauss_values <- dnorm(seq(-gauss_n_indexes_per_side*time_steps,gauss_n_indexes_per_side*time_steps,time_steps), 0, gauss_sd)
  
  # Remove keypresses that are outside relevant time window (time_min - time_max)
  keypresses_to_remove <- time_keypresses < time_min | time_keypresses > time_max
  
  if (sum(keypresses_to_remove) > 0)
  {
    ids <- ids[! keypresses_to_remove]
    time_keypresses <- time_keypresses[! keypresses_to_remove]
    warning(paste0(sum(keypresses_to_remove)," keypresses were removed because they were not within relevant time window (time_min, time_max) == (",time_min,", ",time_max,")"))
  }
  
  # Map each keypress time to an array index
  # ! keep in sync with bootstrap_critical_cutoffs !
  index_keypresses <- round((time_keypresses - time_min) / time_steps)
  
  time_keypresses_rounded <- (index_keypresses * time_steps) + time_min
  
  if (! isTRUE(all.equal(time_keypresses, time_keypresses_rounded)) )
  {
    warning(paste0("Some keypress times were rounded such that they are a multiple of time_min + n * time_steps. Rounded values can be found in variable time_keypresses_rounded. Original values are in variable time_keypresses."))
  }
  
  
  out <- list(
    ids                       =  ids,
    time_keypresses           =  time_keypresses,
    time_keypresses_rounded   =  time_keypresses_rounded,
    time_min                  =  time_min,
    time_max                  =  time_max,
    time_steps                =  time_steps,
    index_time_min            =  0,
    index_time_max            =  index_time_max,
    index_keypresses          =  index_keypresses,    
    indexes_gauss_offset      =  indexes_gauss_offset,
    gauss_n_indexes_per_side  =  gauss_n_indexes_per_side,
    gauss_values              =  gauss_values
  )
  
  class(out) <- "segmag"
  
  # Required at multiple places (e.g., plot, get_eb_times), therefore calculated once directly in constructor
  out$data = calc_segmentation_magnitude(out)
  
  return( out )
}

is.segmag <- function(segmag)
{
  return( class(segmag) == "segmag" )
}

#' Plot segmentation magnitude
#' 
#' Draws a plot depicting the segmentation magnitude resulting from overlaid
#' Gaussians for each key press time across participants. If segmag_substract
#' is defined then the difference in segmentation magnitude of segmag - segmag_substract
#' is drawn.
#' 
#' @param x object of class \code{\link{segmag}}
#' @param cutoffs numeric vector of critical cutoffs that are drawn as horizontal
#'                red lines. Use \code{\link{bootstrap_critical_cutoffs}} in order to
#'                determine the cutoffs for a specific segmag object.
#' @param eb_times numeric vector of event boundary times to highlight in the plot
#' @param segmag_substract object of class \code{\link{segmag}}. If this value is set
#'                         than the difference in segmentation magnitude of
#'                         segmag - segmag_substract is drawn.
#' @param save_as_png string, optional name of file where to save plot (.png is added automatically)
#' @param ... paramters passed to generic plot function
#' @seealso \code{\link{segmag}}
#' @export
plot.segmag <- function(x, cutoffs = NULL, eb_times = NULL, segmag_substract = NULL, save_as_png = NULL, ...)
{
  if (! is.segmag(x)) stop("x must be an object of class segmag")
  
  segmag <- x
  
  cutoffs <- unlist(cutoffs) # convert critical cutoffs returned by bootstrap function when segmag_substract was defined to a vector
  if (! is.null(cutoffs) && ! is.numeric(cutoffs)) stop("cutoffs must be numeric or NULL")
  if (! is.null(eb_times) && ! is.numeric(eb_times)) stop("eb_times must be numeric or NULL")
    
  if (! is.null(save_as_png))
  {
    png(paste(save_as_png,".png",sep=""),960,480)
  }
  
  # Defined here, such that ylab can be replaced in the following if clause
  plot.args <- list(...)
  
  # Calculate difference of segmentation magnitudes if segmag_substract is given
  # and replace values in segmag such that same functions can be used that are
  # used with a normal segmag object without segmag_substract
  if (! is.null(segmag_substract))
  {
    if (! is.segmag(segmag_substract)) stop("segmag_substract must be an object of class segmag")
    if (segmag$time_min != segmag_substract$time_min || segmag$time_max != segmag_substract$time_max || segmag$time_steps != segmag_substract$time_steps) stop("time_min, time_max, and time_steps of segmag and segmag_substract must be equal")
    if (length(segmag$data$time) != length(segmag_substract$data$time) || sum(segmag$data$time != segmag_substract$data$time) > 0) stop("segmag$data$time and segmag_substract$data$time must be equal")
    segmag$data$segmentation_magnitude <- segmag$data$segmentation_magnitude - segmag_substract$data$segmentation_magnitude
    if (is.null(plot.args$ylab)) plot.args$ylab = "Diff in Segmentation Magnitude"
  }

  if (is.null(plot.args$ylim)) plot.args$ylim <- c(min(segmag$data$segmentation_magnitude,cutoffs)-0.1, max(segmag$data$segmentation_magnitude,cutoffs)+0.1)
  if (is.null(plot.args$xlab)) plot.args$xlab = "Time"
  if (is.null(plot.args$ylab)) plot.args$ylab = "Segmentation Magnitude"

  do.call("plot", c(segmag$data$segmentation_magnitude ~ segmag$data$time, plot.args))
  
  for (co in cutoffs)
  {
    abline(h=co, lwd=2, col="red")
  }
  for (eb in eb_times)
  {
    abline(v=eb, lwd=2, col="green")
  }
  
  if (! is.null(save_as_png))
  {
    dev.off()
  }
}