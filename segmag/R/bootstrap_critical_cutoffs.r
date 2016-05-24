#' Bootstrap critical segmentation magnitude
#' 
#' Bootstraps critical segmentation magnitude values for a \code{\link{segmag}}
#' object under the null hypothesis that all key presses were randomly
#' distributed (uniformly) across the experiment (time_min to time_max).
#' 
#' During each bootstrapping iteration, the key presses are randomly distributed
#' (drawn from uniform distribution ranging from time_min to time_max). Then,
#' segmentation magnitude is calculated with those random key press times (note
#' that ids are retained, that is each participant "makes" the same amount of
#' key presses as in the original experiment). The local maxima in segmentation
#' magnitude resulting from the random key press times are ordered according to
#' their size. The largest maximum is kept.
#' 
#' The function returns the \code{critical_probs} quantiles of the vector of those
#' largest maxima obtained across \code{n_bootstap} iterations
#' 
#' This function can also be used to bootstrap the critical maxima and minima
#' cutoffs of a difference function of two segmag objects. To do so, segmag and
#' segmag_substract must be defined. All values will be related to the difference
#' of segmag - segmag_substract (Keypress times in segmag and segmag_substract are
#' randomizes independently).
#' 
#' @param segmag object of class \code{\link{segmag}}
#' @param n_bootstrap numeric, number of bootstrap iterations
#' @param critical_probs numeric vector of probabilities, e.g. c(.95,.99)
#' @param segmag_substract object of class \code{\link{segmag}}. If this value is set
#'                         than keypress times in segmag and segmag_substract are
#'                         both randomized and the critical_cutoffs relate to
#'                         the difference of the segmentation magnitude in segmag
#'                         minus the segmentation magnitude in segmag_substract
#' @param visualize logical, visualize ordered maxima (and minima if segmag_substract
#'                  is set) of bootstrapping iterations (Note: Enabling this option
#'                  might require a lot of RAM with large data sets or large values
#'                  of n_bootstrap)
#' @param save_as character, filename where to save raw bootstrapping data and
#'                plot (optional)
#' @return critical segmentation magnitudes; If segmag_substract is NULL, then
#'         the return value is a numeric vector. Otherwise a list with critical
#'         maxima cutoffs and critical minima cutoffs is returned.
#' @seealso \code{\link{get_eb_times}}, \code{\link{segmag}}
#' @examples #see ?segmag for an example
#' @export
bootstrap_critical_cutoffs <- function(segmag, n_bootstrap, critical_probs, segmag_substract = NULL, visualize = FALSE, save_as = NULL)
{
  if (! is.segmag(segmag)) stop("segmag must be an object of class segmag")
  
  if (! is.null(segmag_substract) && ! is.segmag(segmag_substract)) stop("segmag_substract must be an object of class segmag")
  if (! is.null(segmag_substract))
  {
    if (segmag$time_min != segmag_substract$time_min || segmag$time_max != segmag_substract$time_max || segmag$time_steps != segmag_substract$time_steps) stop("time_min, time_max, and time_steps of segmag and segmag_substract must be equal")
    if (length(segmag$data$time) != length(segmag_substract$data$time) || sum(segmag$data$time != segmag_substract$data$time) > 0) stop("segmag$data$time and segmag_substract$data$time must be equal")
  }
  
  if (! is.numeric(n_bootstrap) || length(n_bootstrap) != 1) stop("n_bootstrap must be numeric and of length 1")
  if (! is.numeric(critical_probs)) stop("critical_probs must be numeric")
  if (! is.logical(visualize) || length(visualize) != 1) stop("visualize must be logical and of length 1")
  if (! is.null(save_as) && (! is.character(save_as) || length(save_as) != 1)) stop("save_as must be a character and of length 1 or NULL")
  
  iteration_ordered_maxima <- list()
  iteration_ordered_minima <- list()
  maxima <- numeric(n_bootstrap)
  minima <- numeric(n_bootstrap)
  
  progress_bar <- txtProgressBar(min = 1, max = n_bootstrap, initial = 1, style = 3)
  for (i in 1:n_bootstrap)
  {
    segmag$time_keypresses <- runif(length(segmag$time_keypresses), segmag$time_min, segmag$time_max)
    
    # ! keep in sync with segmag (class definition) + below !
    segmag$index_keypresses <- round((segmag$time_keypresses - segmag$time_min) / segmag$time_steps)
    
    seg_magnitude <- calc_segmentation_magnitude(segmag)
    
    # If segmag_substract is given, than the segmentation magnitude resulting from
    # random keypress times in segmag_substract is substracted from the segmentation
    # magnitude calculated above
    if (! is.null(segmag_substract))
    {
      segmag_substract$time_keypresses <- runif(length(segmag_substract$time_keypresses), segmag_substract$time_min, segmag_substract$time_max)
      
      # ! keep in sync with segmag (class definition) + above !
      segmag_substract$index_keypresses <- round((segmag_substract$time_keypresses - segmag_substract$time_min) / segmag_substract$time_steps)
      
      seg_magnitude_substract <- calc_segmentation_magnitude(segmag_substract)
      
      if (length(seg_magnitude$time) != length(seg_magnitude_substract$time) || sum(seg_magnitude$time != seg_magnitude_substract$time) > 0) stop("Please contact the package maintainer and report this bug including your data set (id: bootstrap_critical_cutoffs_seg_time). Thanks!")
      seg_magnitude$segmentation_magnitude <- seg_magnitude$segmentation_magnitude - seg_magnitude_substract$segmentation_magnitude
    }
    
    cur_maxima <- seg_magnitude$segmentation_magnitude[flag_maxima_positions(seg_magnitude$segmentation_magnitude)]
    cur_ordered_maxima <- cur_maxima[order(cur_maxima, decreasing=TRUE)]
    if (visualize) iteration_ordered_maxima[[i]] <- cur_ordered_maxima
    maxima[i] <- cur_ordered_maxima[1]
    #print(paste(i, ". max=", maxima[i], sep=""))
    
    if (! is.null(segmag_substract))
    {
      cur_minima <- seg_magnitude$segmentation_magnitude[flag_minima_positions(seg_magnitude$segmentation_magnitude)]
      cur_ordered_minima <- cur_minima[order(cur_minima, decreasing=FALSE)]
      if (visualize) iteration_ordered_minima[[i]] <- cur_ordered_minima
      minima[i] <- cur_ordered_minima[1]
      #print(paste(i, ". min=", minima[i], sep=""))
    }
    
    setTxtProgressBar(progress_bar, i)
  }
  close(progress_bar)
  critical_cutoffs_max <- quantile(maxima, critical_probs)
  critical_cutoffs_min <- NA
  if (! is.null(segmag_substract))
  {
    critical_cutoffs_min <- quantile(minima, 1-critical_probs)
  }
  
  if (visualize)
  {
    # use segmag for difference scores if segmag_substract is defined
    if (! is.null(segmag_substract))
    {
      segmag$data$segmentation_magnitude <- segmag$data$segmentation_magnitude - segmag_substract$data$segmentation_magnitude
    }
    
    orig_maxima <- segmag$data$segmentation_magnitude[flag_maxima_positions(segmag$data$segmentation_magnitude)]
    orig_maxima_ordered <- orig_maxima[order(orig_maxima, decreasing=TRUE)]
    if (! is.null(save_as))
    {
      png(paste(save_as,".png",sep=""), 960, 480)
    }
    plot(NULL, main="ordered maxima", ylab="Segmentation Magnitude", xlim = c(1, length(orig_maxima_ordered)), ylim = c(min(orig_maxima_ordered, min(maxima))-0.1 , max(orig_maxima_ordered, max(maxima))+0.1))
    for (i in 1:n_bootstrap) {
      points(iteration_ordered_maxima[[i]], type="l")
    }
    points(orig_maxima_ordered, type="l", lwd=2, col="blue")
    for (cs in critical_cutoffs_max)
    {
      abline(h=cs, lwd=2, col="red")
    }
    if (! is.null(save_as))
    {
      dev.off()
    }
    
    if (! is.null(segmag_substract))
    {
      orig_minima <- segmag$data$segmentation_magnitude[flag_minima_positions(segmag$data$segmentation_magnitude)]
      orig_minima_ordered <- orig_minima[order(orig_minima, decreasing=FALSE)]
      if (! is.null(save_as))
      {
        png(paste(save_as,"_ordered_minima.png",sep=""), 960, 480)
      }
      plot(NULL, main="ordered minima", ylab="Segmentation Magnitude", xlim = c(1, length(orig_minima_ordered)), ylim = c(min(orig_minima_ordered, min(minima))-0.1 , max(orig_minima_ordered, max(minima))+0.1))
      for (i in 1:n_bootstrap) {
        points(iteration_ordered_minima[[i]], type="l")
      }
      points(orig_minima_ordered, type="l", lwd=2, col="blue")
      for (cs in critical_cutoffs_min)
      {
        abline(h=cs, lwd=2, col="red")
      }
      if (! is.null(save_as))
      {
        dev.off()
      }
    }
  }
  
  if (! is.null(save_as))
  {
    save(list = c("iteration_ordered_maxima","iteration_ordered_minima","critical_cutoffs_max","critical_cutoffs_min"), file= paste(save_as,".RData",sep=""))
  }
  
  if (is.null(segmag_substract))
  {
    return(critical_cutoffs_max)
  }
  else
  {
    return(list(critical_cutoffs_max=critical_cutoffs_max, critical_cutoffs_min=critical_cutoffs_min))
  }
}