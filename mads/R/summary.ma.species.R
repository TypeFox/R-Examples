#' Print a summary of an element of a multi-analysis result correcponding to 
#' a single species included in the analyses.
#' 
#' Provides a summary of the fitted detection probability model
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error. What is printed depends
#' on the corresponding call to summary.
#' 
#' @export summary ma.species
#' @method summary ma.species
#' @aliases summary.ma.species
#' @param object a summary of \code{ma} model object
#' @param species optional character value giving the species name, solely for 
#'   display purposes
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @author Laura Marshall
#' @seealso \code{\link{summary.ma}}
#' @keywords utility
#' @importFrom stats sd
summary.ma.species <- function (object, species=NULL, ...){

  print.tables <- function(object){
    cat("\nBootstrap summary statistics:\n")
    print(object$summary)
    if("N" %in% names(object)){
      cat("\nAbundance:\n")
      print(object$N)
    }
    cat("\nDensity:\n")
    print(object$D)
  }
  
  #Display title line
  if(!is.null(species)){
    cat("\nBootstrap summary for species  : ",species,"\n")
  }else{
    cat("\nBootstrap summary for individual species in a multi-analysis object\n")
  }  
  #Display ddf summary
  cat("\nDetection function model summary\n")
  cat("\nModel Selection:\n")
  print(object$ddf$convergence)
  model.names <- dimnames(object$ddf$convergence)[2][[1]]
  criteria <- names(object$ddf[[model.names[1]]])[2]
  #cat(paste("\nSummary of ", criteria, " values:\n"), sep = "")
  #for(m in seq(along = model.names)){
    #cat(paste("\n", model.names[m], ":\n", sep = ""))
    #print(summary(object$ddf[[model.names[m]]][[criteria]]))    
  #} 
  #cat("\nDetection Function Parameters\n")
  cat("\nModel Summaries\n")
  for(m in seq(along = model.names)){
    cat("\nModel name: ", model.names[m], "\n")
    cat("\nDetection function:\n")
    #print(model.description(get(model.names[m])))
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
  cat("\nDensity / Abundance Summaries\n")
  #Display density (and abundance and expected cluster size tables)
  if(is.null(object$clusters)){
    cat("\nSummary for individuals\n")
    print.tables(object$individuals)
  }else{
    cat("\nSummary for clusters\n")
    print.tables(object$clusters)
    cat("\nSummary for individuals\n")
    print.tables(object$individuals)
    cat("\nExpected cluster size\n\n")
    print(object$Expected.S)
  } 
  invisible()
}