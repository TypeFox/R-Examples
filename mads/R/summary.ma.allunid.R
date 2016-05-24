#' Summary of multi-analysis object
#' 
#' Provides asummary of the fitted detection probability model
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error for all species.
#' 
#' @export summary ma.allunid
#' @method summary ma.allunid
#' @aliases summary.ma.allunid
#' @param object a \code{ma} model object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return list of extracted and summarized objects
#' @note This function is called by the generic function \code{summary} for any
#'   \code{ma} object.  
#' @author Laura Marshall
#' @keywords utility
summary.ma.allunid <- function(object, ...){
  cat("\nUnidentified Model Summaries")
  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  species.name <- names(object)   
  for(sp in seq(along = species.name)){
    summary(object[[sp]], species = species.name[sp])
  }     
  invisible(object)
}