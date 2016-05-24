#' @title Print method for an object of class \code{mbl}
#' @description Prints the content of an object of class \code{mbl}
#' @aliases print.mbl
#' @usage \method{print}{mbl}(x, ...)
#' @param x an object of class \code{mbl} (as returned by the \code{mbl} function). 
#' @param ... arguments to be passed to methods (not yet functional).
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @export

print.mbl <- function(x, ...){
  
  object <- x
  
  if(!is.null(object$nnValStats)){
    nnValStats <- object$nnValStats
  } else {nnValStats <- NULL}
  if(!is.null(object$localCrossValStats)){
    localCrossValStats <- object$localCrossValStats
  } else {localCrossValStats <- NULL}
  if(!is.null(object$YuPredictionStats)){
    YuPredictionStats <- object$YuPredictionStats
  } else {YuPredictionStats <- NULL}
  
  cat("\n")
  cat("Call:", "\n\n")
  print(object$call)
  cat("\n")
  cat("---------------------------------------------------------", "\n")
  cat("\n", "Total number of samples predicted:",object$totalSamplesPredicted,  "\n\n")  
  cat("---------------------------------------------------------", "\n")
  cat("\n", "Total number of PCs used:",object$pcAnalysis$n.componentsUsed,  "\n\n")  
  cat("---------------------------------------------------------", "\n")
  
  if(!is.null(nnValStats)){
    cat("\n", "Nearest neighbor validation statistics", "\n\n")
    print(nnValStats, digits = 3)
    cat("---------------------------------------------------------", "\n")}
  
  if(!is.null(localCrossValStats)){
    cat("\n", "Average statistics of the local leave-group-out","\n", "cross-validation", "\n\n")
    print(localCrossValStats, digits = 3)
    cat("---------------------------------------------------------", "\n")}
  
  if(!is.null(YuPredictionStats)){
    cat("\n", "Statistics of the prediction of Yu", "\n\n")  
    print(YuPredictionStats, digits = 3)}
}
