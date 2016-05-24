#' Checks whether the method is implemented and whether all dependencies are
#' installed
#'
#' @param method Method to use.
#'
checkDependencies <- function(method){
  if(method %in% c("ICP", "hiddenICP")){
    missingDependenciesMessage("InvariantCausalPrediction", method)
  }else if(method %in% c("pc", "LINGAM", "ges", "gies", "rfci")){
    missingDependenciesMessage("pcalg", method)
  }else if(method == "CAM"){
    missingDependenciesMessage("CAM", method)
  }else if(method == "regression"){
    missingDependenciesMessage("glmnet", method)
  }else if(method == "bivariateCAM"){
    missingDependenciesMessage("mgcv", method)
  }else if(method == "bivariateANM"){
    missingDependenciesMessage("kernlab", method)
    missingDependenciesMessage("mgcv", method)
  }else{
    if(method != "backShift")
      stop(paste("Method", method, "not (yet?) implemented."))
  }
}


missingDependenciesMessage <- function(package, method){
  if(!requireNamespace(package, quietly = TRUE)){
    stop(paste("The package '", package, "' is needed for ", method," to 
work. Please install it.", sep=""),
         call. = FALSE)
  }
}