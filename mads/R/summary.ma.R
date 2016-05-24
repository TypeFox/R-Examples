#' Summary of multi-analysis object
#' 
#' Provides asummary of the fitted detection probability model
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error for all species.
#' 
#' @export summary ma
#' @method summary ma
#' @aliases summary.ma
#' @param object an object of class \code{ma} 
#' @param glossary a \code{ma} model object would like a glossary to 
#'   be displayed 
#' @param description boolean if you would like  
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return list of extracted and summarized objects
#' @note This function is called by the generic function \code{summary} for any
#'   \code{ma} object.  
#' @author Laura Marshall
#' @keywords utility
summary.ma <- function(object, description = FALSE, glossary = FALSE,  ...){
  if(description){
    #summary description
    object.description()    
  }
  if(glossary){
    #summary description
    glossary()    
  }
  cat("\nMULTI-ANALYSIS SUMMARY") 
  cat("\n----------------------\n")  
  for(ele in seq(along = object)){
    summary(object[[ele]])
  } 
  invisible(object)
}