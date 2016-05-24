#' Summary CEC
#' @export
#' @rdname summary.cec
#' @method summary Rcpp_CecModel 
#' 
#' @title summary
#' 
#' @description Print detailed information about CEC model object
#'
#' @docType methods
#'
#' @param object CEC model object.
#' @param ... other arguments not used by this method.
#' 
summary.Rcpp_CecModel <- function(object, ...) {
  print(object)
  
  if(isParameterOn(object$iterations)){
    print("Iterations: ")
    print(object$iterations)
  }
  if(isParameterOn(object$logEnergy)){
    print("Energy for every iteration: ")
    print(object$logEnergy)
  }
  if(isParameterOn(object$logNumberOfClusters)){
    print("Number of clusters for every iteration: ")
    print(object$logNumberOfClusters)
  }
}

show.Rcpp_CecModel <- function(object) {
    summary(object)
}

isParameterOn <- function(x) {
  return(length(x) != 0)
}
