#' Check the validity of the stability selection inputs.
#'
#' @param p Number of variables
#' @param q Stability selection parameter
#' @param nodewise Do stability selection nodewise?
#' @param threshold Threshold for stability selection.
#' @param sampleSettings The proportion of unique settings to resample for each 
#'  resample; has to be in [0,1].
#' @param sampleObservations The fraction of all samples to retain when 
#'  subsampling (no replacement); has to be in [0,1].
#'
stabilitySelectionChecks <- function(p, q, nodewise, threshold, sampleSettings,
                                     sampleObservations){
  # check validity of q
  if(q<1){
    qmin <- if(!nodewise) 1/((2*threshold-1)*(p^2-p)) else 1/(2*threshold-1)
    stop(paste("Number of selected edges in each subsample is less than 1. 
               Need to increase EV (expected number of false positives) to 
               at least ", signif(qmin,3), " to have the chance of getting 
               an interesting result", 
               if(nodewise) " (or switch 'nodewise' to FALSE)." else "."),
         sep="")
  }
  # check remaining stability selection inputs
  if(sampleSettings<=0) stop("sampleSettings needs to be positive")
  if(sampleObservations<=0) stop("sampleObservations needs to be positive")
  if(sampleSettings>1) stop("sampleSettings needs to be at most 1")
  if(sampleObservations>1) stop("sampleObservations needs to be at most 1")
}