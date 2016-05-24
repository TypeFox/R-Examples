#' Summary of multi-analysis object
#' 
#' Provides asummary of the fitted detection probability model
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error for all species.
#' 
#' @export summary ma.unid
#' @method summary ma.unid
#' @aliases summary.ma.unid
#' @param object an object of class \code{ma.unid} 
#' @param species optional character value giving the species name, solely for 
#'   display purposes
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return list of extracted and summarized objects
#' @note This function is called by the generic function \code{summary} for any
#'   \code{ma} object.  
#' @author Laura Marshall
#' @keywords utility
#' @importFrom stats sd
summary.ma.unid <- function(object, species = NULL, ...){
  if(!is.null(object)){
    #Display title line
    if(!is.null(species)){
      cat("\nBootstrap summary for code  : ",species,"\n")
    }else{
      cat("\nBootstrap summary for individual unidentified code in a multi-analysis object\n")
    }  
    #Display ddf summary
    cat("\nDetection function model summary\n")
    cat("\nModel Selection:\n")
    print(object$ddf$convergence)
    model.names <- dimnames(object$ddf$convergence)[2][[1]]
    criteria <- names(object$ddf[[model.names[1]]])[2]
    cat("\nModel Summaries\n")
    for(m in seq(along = model.names)){
      cat("\nModel name: ", model.names[m], "\n")
      cat("\nDetection function:\n")
      cat(object$ddf[[model.names[m]]]$model.description)
      cat("\n")
      selected <- FALSE
      if(!is.null(object$ddf[[model.names[m]]]$ds.param)){
        cat("\nParameter estimates (dsmodel):\n")
        if(class(object$ddf[[model.names[m]]]$ds.param) == "matrix"){
          param.estimates <- apply(object$ddf[[model.names[m]]]$ds.param, 2, mean)
          param.se <- apply(object$ddf[[model.names[m]]]$ds.param, 2, sd)
          print(array(c(param.estimates, param.se), dim=c(length(param.estimates),2), dimnames=list(dimnames(object$ddf[[model.names[m]]]$ds.param)[[2]], c("Estimate", "se"))))
        }else if(class(object$ddf[[model.names[m]]]$ds.param) == "numeric"){
          param.estimates <- object$ddf[[model.names[m]]]$ds.param
          param.se <- rep(NA, length(object$ddf[[model.names[m]]]$ds.param))
          print(array(c(param.estimates, param.se), dim=c(length(param.estimates),2), dimnames=list(names(object$ddf[[model.names[m]]]$ds.param), c("Estimate", "se"))))
        }else{
          cat("\nModel never selected\n")
          selected <- FALSE
        }
      }
      if(!is.null(object$ddf[[model.names[m]]]$mr.param)){
        cat("\nParameter estimates (mrmodel):\n")
        if(class(object$ddf[[model.names[m]]]$mr.param) == "matrix"){
          param.estimates <- apply(object$ddf[[model.names[m]]]$mr.param, 2, mean)
          param.se <- apply(object$ddf[[model.names[m]]]$mr.param, 2, sd)
          print(array(c(param.estimates, param.se), dim=c(length(param.estimates),2), dimnames=list(dimnames(object$ddf[[model.names[m]]]$mr.param)[[2]], c("Estimate", "se"))))
        }else if(class(object$ddf[[model.names[m]]]$mr.param) == "numeric"){
          param.estimates <- object$ddf[[model.names[m]]]$mr.param
          param.se <- rep(NA, length(object$ddf[[model.names[m]]]$mr.param))
          print(array(c(param.estimates, param.se), dim=c(length(param.estimates),2), dimnames=list(names(object$ddf[[model.names[m]]]$mr.param), c("Estimate", "se"))))
        }else{
          cat("\nModel never selected\n")
          selected <- FALSE
        }        
      }
      if(selected){
        cat(paste("\nSummary of ", criteria, " values:\n"), sep = "")
        print(summary(object$ddf[[model.names[m]]][[criteria]]))  
      }         
    }
   
  }
  invisible(object)
}