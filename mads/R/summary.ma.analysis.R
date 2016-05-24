#' Summary of multi-analysis object
#' 
#' Provides asummary of the fitted detection probability model
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error for all species.
#' 
#' @export summary ma.analysis
#' @method summary ma.analysis
#' @aliases summary.ma.analysis
#' @param object a \code{ma} model object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return list of extracted and summarized objects
#' @note This function is called by the generic function \code{summary} for any
#'   \code{ma} object.  
#' @author Laura Marshall
#' @keywords utility
summary.ma.analysis <- function(object, ...){
  cat("\nSummary of Analysis Details")
  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("\nBootstrap resample implemented:", object$bootstrap)
  if(!is.null(object$covariate.uncertainty)){
    cat("\nThe following covariates were resampled parametrically and/or had a correction factor applied:", paste(object$covariate.uncertainty$variable.name, collapse = ", ", sep=""))
  }
  if(object$bootstrap | !is.null(object$covariate.uncertainty)){
    cat("\nNumber of resamples:", object$n)
  }
  if(object$unidentified.species){
    cat("\n\nUnidentified species codes were included in these analyses. \nThey were prorated as follows:")
    for(sp in seq(along = object$species.code.definitions)){
      if(length(object$species.code.definitions[[sp]]) > 1){
        cat("\n   Unidentified code ",names(object$species.code.definitions)[sp]," was prorated to species codes ", paste(object$species.code.definitions[[sp]], collapse = ", "), sep="")
      }
    }
  }
  if(length(which(sapply(object$model.name, length) > 1)) > 1){
    cat("\n\nModel uncertainty was included in these analyses. See species results for convergence, selection and model summaries.")
  }
  cat("\n\nData details:")
  cat("\n   Clusters: ", object$clusters, sep="")
  cat("\n   Double observer: ", object$double.observer, "\n", sep="")   
  invisible(object)
}